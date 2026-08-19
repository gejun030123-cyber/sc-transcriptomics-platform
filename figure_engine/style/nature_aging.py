"""Nature Aging 的轻量风格变体。"""

from .nature import NatureStyle


class NatureAgingStyle(NatureStyle):
    name = 'nature_aging'
    signal_blue = '#496D88'
    signal_teal = '#5E8E82'
    signal_red = '#AA625D'
    deg_palette = {'NS': '#DADDE0', 'Up': '#AA625D', 'Down': '#496D88'}
