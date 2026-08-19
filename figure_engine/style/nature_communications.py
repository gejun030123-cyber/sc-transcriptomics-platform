"""Nature Communications 的轻量风格变体。"""

from .nature import NatureStyle


class NatureCommunicationsStyle(NatureStyle):
    name = 'nature_communications'
    signal_blue = '#3F739D'
    signal_red = '#B45F5F'
    deg_palette = {'NS': '#D9DDE2', 'Up': '#B45F5F', 'Down': '#3F739D'}
