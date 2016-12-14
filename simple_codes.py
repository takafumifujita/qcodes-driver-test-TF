# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 11:31:28 2016

@author: tfujita
"""

#%% First load
import qcodes
import qcodes.instrument_drivers.oxford.ILM200 as ILM200

server_name = None

helium = ILM200.OxfordInstruments_ILM200(name='helium', address='ASRL4::INSTR', server_name=server_name)

#%% Reloading

helium.close()

from imp import reload
reload(qcodes.instrument_drivers.oxford.ILM200)
import qcodes.instrument_drivers.oxford.ILM200 as ILM200

server_name = None

helium = ILM200.OxfordInstruments_ILM200(name='helium', address='ASRL4::INSTR', server_name=server_name)
