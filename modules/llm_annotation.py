"""De-identified cluster-level LLM annotation adapter.

This adapter deliberately sends only a compact, cluster-level marker summary
to the configured model provider.  It never receives AnnData, expression
values, cell barcodes, donor/sample metadata, project paths, or free-form
annotation comments.
"""
from __future__ import annotations

import hashlib
import json
import re
from typing import Any, Dict, Iterable, List, Optional

from modules.ai_config import get_effective_ai_config


LLM_ANNOTATION_PROMPT_VERSION = "cluster-marker-json-v1"
MAX_CONTEXT_LENGTH = 500
MAX_CLUSTERS = 40
MAX_MARKERS_PER_CLUSTER = 12
# 每个 cluster 需要输出 cell_type + 240 字符 rationale + 240 字符 review_note
# 以及 JSON 结构字段。中文按“1 字符 ≈ 1 token”保守估算，避免 30 簇的批次
# 在 4096 token 处被静默截断（截断会被误判为“模型没按契约返回”）。
TOKENS_PER_CLUSTER = 320
MIN_RESPONSE_TOKENS = 2048
MAX_RESPONSE_TOKENS = 16000


def _response_token_budget(n_clusters: int) -> int:
    """Return a response budget that can actually cover the requested batch."""
    try:
        count = max(1, int(n_clusters))
    except (TypeError, ValueError):
        count = 1
    return int(min(MAX_RESPONSE_TOKENS, max(MIN_RESPONSE_TOKENS, count * TOKENS_PER_CLUSTER)))


class LLMAnnotationError(RuntimeError):
    """Raised when the configured LLM cannot return a valid annotation."""


class LLMResponseError(LLMAnnotationError):
    """Provider answered, but the reply did not satisfy the JSON contract.

    Distinguished from transport/configuration failures so a truncated or
    incomplete response can be retried with a smaller batch instead of
    failing the whole annotation task.
    """


def _provider(settings: Dict[str, str]) -> str:
    provider = str(settings.get("provider") or "auto").strip().lower()
    if provider == "anthropic":
        return "anthropic"
    if provider == "openai":
        return "openai"
    url = str(settings.get("api_url") or "").lower()
    return "anthropic" if "anthropic" in url or "claude" in url else "openai"


def _anthropic_endpoint(base_url: str) -> str:
    base_url = str(base_url or "").rstrip("/")
    if base_url.endswith("/v1/messages") or base_url.endswith("/messages"):
        return base_url
    if base_url.endswith("/v1"):
        return base_url + "/messages"
    return base_url + "/v1/messages"


def _openai_endpoint(base_url: str) -> str:
    base_url = str(base_url or "").rstrip("/")
    if base_url.endswith("/chat/completions"):
        return base_url
    if base_url.endswith("/v1"):
        return base_url + "/chat/completions"
    return base_url + "/v1/chat/completions"


def _plain_text(value: Any, *, limit: int) -> str:
    text = " ".join(str(value or "").split())
    return text[:limit]


def _generic_tissue_context(value: Any) -> str:
    """Keep optional biological context small and reject obvious identifiers."""
    text = _plain_text(value, limit=MAX_CONTEXT_LENGTH)
    if not text:
        return ""
    if re.search(r"[\w.+-]+@[\w.-]+\.[A-Za-z]{2,}", text):
        raise LLMAnnotationError("LLM 组织背景不能包含邮箱或直接联系方式。")
    if re.search(r"(?:^|\s)(?:/home/|/data/|[A-Za-z]:[\\/]|~[\\/])", text):
        raise LLMAnnotationError("LLM 组织背景不能包含服务器路径。")
    return text


def _normalise_cluster_summaries(
    cluster_summaries: Iterable[Dict[str, Any]], *,
    max_clusters: Optional[int] = None,
) -> List[Dict[str, Any]]:
    """Normalise cluster summaries without leaking extra fields.

    ``max_clusters`` is the per-request batch cap used by the payload builder.
    ``None`` (used by the chunked runner) means "normalise everything"; the
    caller is responsible for splitting into per-request batches.
    """
    if max_clusters is not None:
        try:
            max_clusters = int(max_clusters)
        except (TypeError, ValueError):
            max_clusters = 30
        max_clusters = max(1, min(max_clusters, MAX_CLUSTERS))

    normalised: List[Dict[str, Any]] = []
    seen = set()
    for item in cluster_summaries or []:
        if not isinstance(item, dict):
            continue
        cluster = _plain_text(item.get("cluster"), limit=64)
        if not cluster or cluster in seen:
            continue
        seen.add(cluster)
        try:
            n_cells = max(0, int(item.get("n_cells", 0)))
        except (TypeError, ValueError):
            n_cells = 0
        markers = []
        for gene in item.get("top_markers", []) or []:
            gene = _plain_text(gene, limit=64)
            if gene and gene not in markers:
                markers.append(gene)
            if len(markers) >= MAX_MARKERS_PER_CLUSTER:
                break
        normalised.append({
            "cluster": cluster,
            "n_cells": n_cells,
            "top_markers": markers,
            # This is a marker-rule candidate, not a patient or cell-level
            # attribute.  The model is explicitly told to review it critically.
            "marker_panel_candidate": _plain_text(
                item.get("marker_panel_candidate"), limit=120,
            ) or "Unknown",
            "marker_panel_support": _plain_text(
                item.get("marker_panel_support"), limit=300,
            ),
        })
    if not normalised:
        raise LLMAnnotationError("没有可供 LLM 注释的 cluster marker 摘要。")
    if max_clusters is not None and len(normalised) > max_clusters:
        raise LLMAnnotationError(
            f"当前有 {len(normalised)} 个 cluster，超过单批 LLM 注释上限 {max_clusters}；"
            f"请提高 llm_max_clusters（单批最大 {MAX_CLUSTERS}）或使用分批注释。"
        )
    return normalised


def build_llm_annotation_payload(
    cluster_summaries: Iterable[Dict[str, Any]],
    *,
    species: Any = "human",
    tissue_context: Any = "",
    max_clusters: int = 30,
) -> Dict[str, Any]:
    """Build the exact de-identified request payload used for annotation."""
    species = _plain_text(species or "unknown", limit=40).lower() or "unknown"
    if species not in {"human", "mouse", "other", "unknown"}:
        species = "unknown"
    return {
        "prompt_version": LLM_ANNOTATION_PROMPT_VERSION,
        "species": species,
        "tissue_context": _generic_tissue_context(tissue_context),
        "clusters": _normalise_cluster_summaries(
            cluster_summaries, max_clusters=max_clusters,
        ),
    }


SYSTEM_PROMPT = """You are a research-only single-cell RNA-seq cell-type annotation assistant.
You receive a de-identified cluster-level summary. It contains only cluster IDs,
cell counts, short lists of top marker genes, and a preliminary marker-panel
candidate. Do not infer donor identity, clinical status, or a diagnosis. Do not
invent genes that are absent from the supplied summaries.

Return JSON only, with this exact top-level structure:
{"annotations":[{"cluster":"...","cell_type":"...","confidence":"high|medium|low|unknown","rationale":"...","review_note":"..."}],"global_note":"..."}

Return exactly one item for every supplied cluster ID. Prefer the most
specific cell type that any supplied marker support points to, even when the
call is uncertain: mark such calls low/unknown confidence and explain the
ambiguity in review_note. Reserve "Unknown" for clusters where no supplied
marker supports any candidate at all. Keep each rationale under 240
characters, cite only supplied marker genes, and state important ambiguity in
review_note. These are research annotations requiring human review, never
clinical conclusions."""


def _extract_json_object(text: Any) -> Dict[str, Any]:
    text = str(text or "").strip()
    if text.startswith("```"):
        text = re.sub(r"^```(?:json)?\s*", "", text, flags=re.IGNORECASE)
        text = re.sub(r"\s*```$", "", text).strip()
    decoder = json.JSONDecoder()
    for index, char in enumerate(text):
        if char != "{":
            continue
        try:
            value, _ = decoder.raw_decode(text[index:])
        except json.JSONDecodeError:
            continue
        if isinstance(value, dict):
            return value
    raise LLMResponseError("LLM 未返回可解析的 JSON 注释结果（可能被截断）。")


def _validate_llm_response(response: Dict[str, Any], expected_clusters: Iterable[str]) -> Dict[str, Any]:
    expected = [str(cluster) for cluster in expected_clusters]
    expected_set = set(expected)
    annotations = response.get("annotations") if isinstance(response, dict) else None
    if not isinstance(annotations, list):
        raise LLMResponseError("LLM 返回缺少 annotations 列表。")

    by_cluster: Dict[str, Dict[str, str]] = {}
    valid_confidence = {"high", "medium", "low", "unknown"}
    for item in annotations:
        if not isinstance(item, dict):
            continue
        cluster = _plain_text(item.get("cluster"), limit=64)
        if cluster not in expected_set or cluster in by_cluster:
            continue
        label = _plain_text(item.get("cell_type"), limit=120)
        if not label:
            raise LLMResponseError(f"LLM 未为 cluster {cluster} 提供 cell_type。")
        confidence = _plain_text(item.get("confidence"), limit=20).lower()
        if confidence not in valid_confidence:
            confidence = "unknown"
        by_cluster[cluster] = {
            "cluster": cluster,
            "cell_type": label,
            "confidence": confidence,
            "rationale": _plain_text(item.get("rationale"), limit=240),
            "review_note": _plain_text(item.get("review_note"), limit=240),
        }
    missing = [cluster for cluster in expected if cluster not in by_cluster]
    if missing:
        raise LLMResponseError(
            "LLM 返回未覆盖全部 cluster（响应可能被截断）：" + ", ".join(missing[:8])
        )
    return {
        "annotations": [by_cluster[cluster] for cluster in expected],
        "global_note": _plain_text(response.get("global_note"), limit=500),
    }


def _response_text(response: Any, provider: str) -> str:
    try:
        body = response.json()
    except (TypeError, ValueError) as exc:
        raise LLMResponseError("LLM 服务返回了非 JSON 响应（响应体可能被截断）。") from exc
    if provider == "anthropic":
        content = body.get("content", []) if isinstance(body, dict) else []
        text = "".join(
            str(block.get("text", ""))
            for block in content if isinstance(block, dict) and block.get("type") == "text"
        )
    else:
        choices = body.get("choices", []) if isinstance(body, dict) else []
        message = choices[0].get("message", {}) if choices and isinstance(choices[0], dict) else {}
        content = message.get("content", "") if isinstance(message, dict) else ""
        # OpenAI 兼容服务可能返回内容块列表（[{type: text, text: ...}]）。
        if isinstance(content, list):
            text = "".join(
                str(block.get("text", ""))
                for block in content if isinstance(block, dict) and block.get("type") == "text"
            )
        else:
            text = str(content)
    if not text:
        raise LLMResponseError("LLM 服务未返回注释文本（响应可能被截断）。")
    return str(text)


def run_llm_cluster_annotation(
    cluster_summaries: Iterable[Dict[str, Any]],
    *,
    species: Any = "human",
    tissue_context: Any = "",
    max_clusters: int = 30,
) -> Dict[str, Any]:
    """Call the configured provider with cluster-only evidence and validate it.

    High-resolution clustering can produce far more clusters than one request
    can carry.  max_clusters is therefore the per-request batch size: when
    more clusters are supplied, they are annotated in multiple independent
    batches so every cluster still receives a label.  Each batch is a separate
    provider call with its own request hash; a single-batch run keeps the
    legacy string request_hash, multi-batch runs return the list.
    """
    import requests

    settings = get_effective_ai_config()
    api_key = str(settings.get("api_key") or "").strip()
    api_url = str(settings.get("api_url") or "").strip()
    model = str(settings.get("model") or "").strip()
    if not api_key:
        raise LLMAnnotationError("AI API Key 未配置；请先在 AI 设置中保存可用的模型服务。")
    if not api_url or not model:
        raise LLMAnnotationError("AI API 地址或模型名未配置；请先完成 AI 设置。")

    try:
        per_batch = max(1, min(int(max_clusters), MAX_CLUSTERS))
    except (TypeError, ValueError):
        per_batch = 30
    normalised = _normalise_cluster_summaries(cluster_summaries, max_clusters=None)
    batches = [
        normalised[index:index + per_batch]
        for index in range(0, len(normalised), per_batch)
    ]
    provider = _provider(settings)

    def call_batch(batch: List[Dict[str, Any]]) -> Dict[str, Any]:
        payload = build_llm_annotation_payload(
            batch,
            species=species,
            tissue_context=tissue_context,
            max_clusters=per_batch,
        )
        prompt = (
            "Annotate the following de-identified cluster summary. Return JSON only.\n"
            + json.dumps(payload, ensure_ascii=False, separators=(",", ":"))
        )
        # 预算随实际批大小伸缩；固定 4096 token 无法容纳 30 簇的
        # rationale + review_note，截断后会被误判为模型违约。
        response_token_budget = _response_token_budget(len(batch))
        try:
            if provider == "anthropic":
                response = requests.post(
                    _anthropic_endpoint(api_url),
                    headers={
                        "Content-Type": "application/json",
                        "x-api-key": api_key,
                        "Authorization": f"Bearer {api_key}",
                        "anthropic-version": "2023-06-01",
                    },
                    json={
                        "model": model,
                        "system": SYSTEM_PROMPT,
                        "messages": [{"role": "user", "content": prompt}],
                        # 与 OpenAI 分支一致：固定采样温度，保证同一请求可复现。
                        "temperature": 0,
                        "max_tokens": response_token_budget,
                    },
                    timeout=60.0,
                )
            else:
                response = requests.post(
                    _openai_endpoint(api_url),
                    headers={
                        "Content-Type": "application/json",
                        "Authorization": f"Bearer {api_key}",
                    },
                    json={
                        "model": model,
                        "messages": [
                            {"role": "system", "content": SYSTEM_PROMPT},
                            {"role": "user", "content": prompt},
                        ],
                        "temperature": 0,
                        "max_tokens": response_token_budget,
                    },
                    timeout=60.0,
                )
        except requests.RequestException as exc:
            # DNS/连接/超时等传输错误统一转成 LLMAnnotationError，保证
            # annotation.py 的用户友好错误路径生效，而不是暴露原始堆栈。
            raise LLMAnnotationError(f"LLM 服务请求失败（网络/传输错误）：{exc}") from exc
        if getattr(response, "status_code", 500) >= 400:
            raise LLMAnnotationError(
                f"LLM 服务请求失败（HTTP {getattr(response, 'status_code', 'unknown')}）。"
            )
        parsed = _extract_json_object(_response_text(response, provider))
        return _validate_llm_response(
            parsed, [item["cluster"] for item in payload["clusters"]],
        )

    split_events: List[Dict[str, Any]] = []

    def annotate_batch(batch: List[Dict[str, Any]]):
        """Annotate one batch, halving it when the reply breaks the contract.

        Returns ``(result, batches_sent)`` so the recorded request hashes still
        describe exactly what was sent.  Only response-contract failures are
        retried; transport and configuration errors propagate unchanged.
        """
        try:
            return call_batch(batch), [batch]
        except LLMResponseError as exc:
            if len(batch) <= 1:
                raise LLMAnnotationError(
                    "LLM 响应无法满足注释契约（cluster "
                    f"{batch[0].get('cluster') if batch else '?'}）：{exc}"
                ) from exc
            mid = max(1, len(batch) // 2)
            split_events.append({
                "n_clusters": len(batch),
                "split_into": [mid, len(batch) - mid],
                "reason": str(exc),
            })
            left, left_sent = annotate_batch(batch[:mid])
            right, right_sent = annotate_batch(batch[mid:])
            merged = {
                "annotations": list(left["annotations"]) + list(right["annotations"]),
                "global_note": " | ".join(
                    note for note in (left.get("global_note"), right.get("global_note"))
                    if note
                ),
            }
            return merged, left_sent + right_sent

    annotations: List[Dict[str, Any]] = []
    global_notes: List[str] = []
    request_hashes: List[str] = []
    for batch in batches:
        batch_result, sent_batches = annotate_batch(batch)
        annotations.extend(batch_result["annotations"])
        if batch_result.get("global_note"):
            global_notes.append(batch_result["global_note"])
        for sent_batch in sent_batches:
            batch_payload = build_llm_annotation_payload(
                sent_batch, species=species, tissue_context=tissue_context,
                max_clusters=per_batch,
            )
            request_hashes.append(hashlib.sha256(
                json.dumps(batch_payload, ensure_ascii=False, sort_keys=True).encode("utf-8")
            ).hexdigest())

    tissue_context_provided = bool(_generic_tissue_context(tissue_context))
    return {
        "enabled": True,
        "provider": provider,
        "model": model,
        "prompt_version": LLM_ANNOTATION_PROMPT_VERSION,
        "request_hash": (
            request_hashes[0] if len(request_hashes) == 1 else request_hashes
        ),
        "cluster_input": normalised,
        "annotations": annotations,
        "global_note": " | ".join(global_notes),
        "n_batches": len(batches),
        "n_requests_sent": len(request_hashes),
        "batch_splits": split_events,
        "max_clusters_per_batch": per_batch,
        "response_token_budget": _response_token_budget(per_batch),
        "tissue_context_provided": tissue_context_provided,
        "sent_data_policy": (
            "仅发送 cluster ID、细胞数、Top marker 与 marker 规则候选；"
            "不发送表达矩阵、细胞条形码、donor/sample 元数据、绝对路径或注释备注。"
        ),
    }
