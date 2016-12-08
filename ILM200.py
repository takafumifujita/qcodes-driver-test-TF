# OxfordInstruments_ILM200.py class, to perform the communication between the Wrapper and the device
# Takafumi Fujita <t.fujita@tudelft.nl>, 2016
# Guenevere Prawiroatmodjo <guen@vvtp.tudelft.nl>, 2009
# Pieter de Groot <pieterdegroot@gmail.com>, 2009
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

## TODO: remove when finished
# from instrument import Instrument
from time import time, sleep
import visa
# import types
import logging
# import numpy
# import struct


from qcodes import VisaInstrument
from qcodes import validators as vals

class OxfordInstruments_ILM200(VisaInstrument):
    """
    This is the qcodes driver for the Oxford Instruments ILM 200 Helium Level Meter.

    Usage:
    Initialize with
    <name> = instruments.create('name', 'OxfordInstruments_ILM200', address='<Instrument address>')
    <Instrument address> = ASRL4::INSTR
    
    Note: Since the ISOBUS allows for several instruments to be managed in parallel, the command
    which is sent to the device starts with '@n', where n is the ISOBUS instrument number.
    
    """
    def __init__(self, name, address, number=1, **kwargs):
        """
        Initializes the Oxford Instruments ILM 200 Helium Level Meter.

        Input:
            name (string)    : name of the instrument
            address (string) : instrument address
            number (int)     : ISOBUS instrument number
							   (number=1 is specific to the ILM in F008)

        Output:
            None
        """
        logging.debug(__name__ + ' : Initializing instrument')
        super().__init__(name, address, **kwargs)

		# TODO: make write and read command based on IVVI
        self.visa_handle.set_visa_attribute(visa.constants.VI_ATTR_ASRL_STOP_BITS,
											visa.constants.VI_ASRL_STOP_TWO)
        self._address = address
        self._number = number
        self._values = {}
        
        self.add_parameter('level',
						   label='level',
						   get_cmd=self._do_get_level,
						   delay=0.02,
						   units='%')
        # self.add_parameter('status', type=types.StringType,
            # flags=Instrument.FLAG_GET)

        # # Add functions
        # self.add_function('get_all')
        # self.get_all()

	def get_idn(self):
		"""
		"""
		return self._get_version

	def get_all(self):
        """
        Reads all implemented parameters from the instrument,
        and updates the wrapper.

        Input:
            None

        Output:
            None
        """
        logging.info(__name__ + ' : reading all settings from instrument')
        self.get_level()
        self.get_status()

    def _execute(self, message):
        """
        Write a command to the device and read answer.
        
        Input:
            message (str) : write command for the device
        
        Output:
            None
        """
        logging.info(__name__ + ' : Send the following command to the device: %s' %message)
        self.visa_handle.write_raw('@%s%s' %(self._number, message))
        sleep(20e-3) # wait for the device to be able to respond
        result = self._read()
        if result.find('?') >= 0:
            print("Error: Command %s not recognized" %message)
        else:
            return result

	def _read(self):
        # because protocol has no termination chars the read reads the number
        # of bytes in the buffer
        bytes_in_buffer = self.visa_handle.bytes_in_buffer
        # a workaround for a timeout error in the pyvsia read_raw() function
        with(self.visa_handle.ignore_warning(visa.constants.VI_SUCCESS_MAX_CNT)):
            mes = self.visa_handle.visalib.read(
                self.visa_handle.session, bytes_in_buffer)
        mes = mes[0]  # cannot be done on same line for some reason
        # if mes[1] != 0:
        #     # see protocol descriptor for error codes
        #     raise Exception('IVVI rack exception "%s"' % mes[1])
		return mes

	# Functions: Monitor commands
    def _get_version(self):
        """
        Identify the device
        
        Input:
            None
            
        Output:
            None
        """
        logging.info(__name__ + ' : Identify the device')
        return self._execute('V')

	def _do_get_level(self):
        """
        Get Helium level of channel 1.
        Input:
            None
        
        Output:
            result (float) : Helium level
        """
        logging.info(__name__ + ' : Read level of channel 1')
        result = self._execute('R1')
        return float(result.replace('R',''))

    def _do_get_status(self):
        """
        Get status of the device.
        Input:
            None
            
        Output:
            None
        """
        logging.info(__name__ + ' : Get status of the device.')
        result = self._execute('X')
        usage = {
        0 : "Channel not in use",
        1 : "Channel used for Nitrogen level",
        2 : "Channel used for Helium Level (Normal pulsed operation)",
        3 : "Channel used for Helium Level (Continuous measurement)",
        9 : "Error on channel (Usually means probe unplugged)"
        }
        # current_flowing = {
        # 0 : "Curent not flowing in Helium Probe Wire",
        # 1 : "Curent not flowing in Helium Probe Wire"
        # }
        # rate = {
        # 10 : "Helium Probe in FAST rate",
        # 01 : "Helium Probe in SLOW rate"
        # }
        # auto_fill_status = {
        # 00 : "End Fill (Level > FULL)",
        # 01 : "Not Filling (Level < FULL, Level > FILL)",
        # 10 : "Filling (Level < FULL, Level > FILL)",
        # 11 : "Start Filling (Level < FILL)"
        # }
        return usage.get(int(result[1]), "Unknown")

    def remote(self):
        """
        Set control to remote & locked
        
        Input:
            None
            
        Output:
            None
        """
        logging.info(__name__ + ' : Set control to remote & locked')
        self.set_remote_status(1)

    def local(self):
        """
        Set control to local & locked
        
        Input:
            None
            
        Output:
            None
        """
        logging.info(__name__ + ' : Set control to local & locked')
        self.set_remote_status(0)

    def set_remote_status(self, mode):
        """
        Set remote control status.
        
        Input:
            mode(int) :
            0 : "Local and locked",
            1 : "Remote and locked",
            2 : "Local and unlocked",
            3 : "Remote and unlocked",
        
        Output:
            None
        """
        status = {
        0 : "Local and locked",
        1 : "Remote and locked",
        2 : "Local and unlocked",
        3 : "Remote and unlocked",
        }
        logging.info(__name__ + ' : Setting remote control status to %s' %status.get(mode,"Unknown"))
        self._execute('C%s' %mode)

    # Functions: Control commands (only recognised when in REMOTE control)
	def set_to_slow(self):
		"""
		Set helium meter channel 1 to slow mode.
		"""
		logging.info(__name__ + ' : Setting Helium Probe in SLOW rate')
		self._execute('S1')

	def set_to_fast(self):
		"""
		Set helium meter channel 1 to fast mode.
		"""
		logging.info(__name__ + ' : Setting Helium Probe in FAST rate')
		self._execute('T1')