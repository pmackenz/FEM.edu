from .SectionMaterial import *

class ElasticSection(SectionMaterial):
    
    def __init__(self, params={}):
        super(ElasticSection, self).__init__(params)
