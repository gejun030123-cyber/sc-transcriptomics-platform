"""Deterministic colour-vision simulation and palette separation audit."""

from pathlib import Path

import numpy as np


_CVD_MATRICES = {
    'protanopia': np.array([
        [0.152286, 1.052583, -0.204868],
        [0.114503, 0.786281, 0.099216],
        [-0.003882, -0.048116, 1.051998],
    ]),
    'deuteranopia': np.array([
        [0.367322, 0.860646, -0.227968],
        [0.280085, 0.672501, 0.047413],
        [-0.011820, 0.042940, 0.968881],
    ]),
    'tritanopia': np.array([
        [1.255528, -0.076749, -0.178779],
        [-0.078411, 0.930809, 0.147602],
        [0.004733, 0.691367, 0.303900],
    ]),
}


def simulate_rgb(rgb, deficiency='deuteranopia'):
    """Simulate complete dichromacy for an sRGB array in the 0–1 range."""

    try:
        matrix = _CVD_MATRICES[str(deficiency).lower()]
    except KeyError as exc:
        raise ValueError(f'未知色觉模拟类型: {deficiency}') from exc
    values = np.asarray(rgb, dtype=float)
    transformed = values @ matrix.T
    return np.clip(transformed, 0.0, 1.0)


def simulate_image(source, destination, deficiency='deuteranopia'):
    """Write a colour-vision preview for a PNG/TIFF-compatible image."""

    from PIL import Image

    with Image.open(source).convert('RGBA') as image:
        values = np.asarray(image, dtype=float) / 255.0
        values[..., :3] = simulate_rgb(values[..., :3], deficiency)
        output = Image.fromarray(np.rint(values * 255).astype(np.uint8))
        destination = Path(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        output.save(destination)
    return str(destination)


def audit_palette(colors, minimum_distance=18.0):
    """Return worst-case pair separation under three dichromacy simulations."""

    from matplotlib.colors import to_rgb

    rgb = np.asarray([to_rgb(color) for color in colors], dtype=float)
    if len(rgb) < 2:
        return {'pass': True, 'minimum_distance': None, 'deficiency': None, 'pair': None}
    worst = {'minimum_distance': float('inf'), 'deficiency': None, 'pair': None}
    for deficiency in _CVD_MATRICES:
        simulated = simulate_rgb(rgb, deficiency)
        for first in range(len(simulated)):
            for second in range(first + 1, len(simulated)):
                distance = float(np.linalg.norm(simulated[first] - simulated[second]) * 255.0)
                if distance < worst['minimum_distance']:
                    worst.update(
                        minimum_distance=distance,
                        deficiency=deficiency,
                        pair=[first, second],
                    )
    worst['minimum_distance'] = round(worst['minimum_distance'], 2)
    worst['pass'] = bool(worst['minimum_distance'] >= float(minimum_distance))
    return worst
