"""Small deterministic visual signatures for regression tests."""

import hashlib

import numpy as np


def figure_signature(fig, hash_size=16):
    """Return a perceptual hash plus structural facts for a Matplotlib figure."""

    from PIL import Image

    fig.canvas.draw()
    rgba = np.asarray(fig.canvas.buffer_rgba(), dtype=np.uint8)
    image = Image.fromarray(rgba).convert('L')
    resized = image.resize((hash_size + 1, hash_size), Image.Resampling.LANCZOS)
    pixels = np.asarray(resized, dtype=np.int16)
    bits = pixels[:, 1:] > pixels[:, :-1]
    packed = np.packbits(bits.reshape(-1))
    return {
        'dhash': packed.tobytes().hex(),
        'pixel_sha256': hashlib.sha256(rgba.tobytes()).hexdigest(),
        'canvas_px': [int(rgba.shape[1]), int(rgba.shape[0])],
        'axes': len(fig.axes),
        'texts': sum(1 for text in fig.findobj()
                     if hasattr(text, 'get_text') and str(text.get_text()).strip()),
        'collections': sum(len(axis.collections) for axis in fig.axes),
        'lines': sum(len(axis.lines) for axis in fig.axes),
    }


def hash_distance(first, second):
    """Hamming distance between two hexadecimal perceptual hashes."""

    left = bytes.fromhex(first)
    right = bytes.fromhex(second)
    if len(left) != len(right):
        raise ValueError('Visual hash 长度不一致')
    return sum((a ^ b).bit_count() for a, b in zip(left, right))


def compare_signatures(actual, expected, maximum_distance=8):
    distance = hash_distance(actual['dhash'], expected['dhash'])
    structural = all(actual[key] == expected[key] for key in (
        'canvas_px', 'axes', 'texts', 'collections', 'lines'
    ))
    return {
        'pass': bool(structural and distance <= int(maximum_distance)),
        'hash_distance': distance,
        'maximum_distance': int(maximum_distance),
        'structural_match': structural,
    }
