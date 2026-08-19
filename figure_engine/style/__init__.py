"""Nature Portfolio 期刊风格注册表。"""

from .nature import NatureStyle, StyleProfile
from .nature_aging import NatureAgingStyle
from .nature_communications import NatureCommunicationsStyle


def get_style(name='nature'):
    styles = {
        'nature': NatureStyle,
        'nature_communications': NatureCommunicationsStyle,
        'nature_aging': NatureAgingStyle,
    }
    try:
        return styles[str(name).lower()]()
    except KeyError as exc:
        raise ValueError(f'未知 Nature style: {name}') from exc


__all__ = [
    'NatureStyle', 'NatureAgingStyle', 'NatureCommunicationsStyle',
    'StyleProfile', 'get_style',
]
