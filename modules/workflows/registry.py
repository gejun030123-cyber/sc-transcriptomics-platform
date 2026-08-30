"""Declarative WES workflow registry.

This registry describes supported contracts.  It intentionally does not run
shell commands; launching a workflow belongs to the separately configured
Nextflow executor.
"""

from typing import Dict, Optional

from .contracts import WorkflowSpec


WES_WORKFLOW_REGISTRY: Dict[str, WorkflowSpec] = {
    "wes_germline": WorkflowSpec(
        key="wes_germline",
        display_name="WES 胚系 SNV/InDel",
        description="HaplotypeCaller gVCF workflow for single samples or small families.",
        modes=("single", "pedigree"),
        input_types=("fastq", "bam", "cram"),
        status="validation",
    ),
    "wes_somatic": WorkflowSpec(
        key="wes_somatic",
        display_name="WES 肿瘤体细胞 SNV/InDel",
        description="Mutect2 workflow for tumor-normal and restricted tumor-only analyses.",
        modes=("tumor_normal", "tumor_only"),
        input_types=("fastq", "bam", "cram"),
        status="validation",
    ),
    "wes_annotate_only": WorkflowSpec(
        key="wes_annotate_only",
        display_name="WES VCF 标准化与注释",
        description="Validate, normalize and annotate an existing compressed VCF offline.",
        modes=("germline", "somatic"),
        input_types=("vcf",),
        status="validation",
    ),
}

# Future assays can share the same API without changing the WES registry.
WORKFLOW_REGISTRY: Dict[str, WorkflowSpec] = dict(WES_WORKFLOW_REGISTRY)


def get_workflow(key: str) -> Optional[WorkflowSpec]:
    return WORKFLOW_REGISTRY.get(str(key or "").strip())


def list_workflows(assay_type: Optional[str] = None):
    specs = WORKFLOW_REGISTRY.values()
    if assay_type:
        specs = [spec for spec in specs if spec.assay_type == assay_type]
    return [spec.to_dict() for spec in specs]
