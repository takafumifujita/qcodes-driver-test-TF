# OxfordInstruments_IPS120.py class, to perform the communication between the Wrapper and the device
# Takafumi Fujita <t.fujita@tudelft.nl>, 2016
# Mohammad Shafiei <m.shafiei@tudelft.nl>, 2011
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

import logging

from qcodes import VisaInstrument
# from qcodes import validators as vals
from time import time, sleep
import visa
import types
import logging
import numpy
import struct


class OxfordInstruments_IPS120(VisaInstrument):
    """This is the python driver for the Oxford Instruments IPS 120 Magnet Power Supply

    Usage:
    Initialize with
    <name> = instruments.create('name', 'OxfordInstruments_IPS120', address='<Instrument address>')
    <Instrument address> = ASRL1::INSTR

    Note: Since the ISOBUS allows for several instruments to be managed in parallel, the command
    which is sent to the device starts with '@n', where n is the ISOBUS instrument number.

    """

# TODO(copied from original): get doesn't always update the wrapper! (e.g.
# when input is an int and output is a string)

    def __init__(self, name, address, number=2, verbose=1, reset=False, **kwargs):
        """Initializes the Oxford Instruments IPS 120 Magnet Power Supply.

        Input:
            name (string)    : name of the instrument
            address (string) : instrument address
            number (int)     : ISOBUS instrument number

        Output:
            None
        """
        self.verbose = verbose
        logging.debug(__name__ + ' : Initializing instrument')
        super().__init__(name, address, **kwargs)

        self._address = address
        self._number = number
        self._visainstrument = visa.SerialInstrument(self._address)
        self._values = {}
        self._visainstrument.stop_bits = 2

        # Add parameters
        self.add_parameter('mode2', type=types.IntType,
                           flags=Instrument.FLAG_GET,
                           format_map={
                               0: "At rest",
                               1: "Sweeping",
                               2: "Sweep limiting",
                               3: "Sweeping and sweep limiting"})
        self.add_parameter(
            'activity',
            type=types.IntType,
            flags=Instrument.FLAG_GETSET | Instrument.FLAG_GET_AFTER_SET,
            format_map={
                0: "Hold",
                1: "To set point",
                2: "To zero",
                4: "Output clamped"})
        self.add_parameter(
            'switch_heater',
            type=types.IntType,
            flags=Instrument.FLAG_GETSET | Instrument.FLAG_GET_AFTER_SET,
            format_map={
                0: "Off magnet at zero (switch closed)",
                1: "On (switch open)",
                2: "Off magnet at field (switch closed)",
                5: "Heater fault (heater is on but current is low)",
                8: "No switch fitted"})
        self.add_parameter(
            'field_setpoint',
            type=types.FloatType,
            flags=Instrument.FLAG_GETSET | Instrument.FLAG_GET_AFTER_SET,
            minval=-8,
            maxval=8)
        self.add_parameter(
            'sweeprate_field',
            type=types.FloatType,
            flags=Instrument.FLAG_GETSET | Instrument.FLAG_GET_AFTER_SET,
            minval=0,
            maxval=0.524)

        # Find the F field limits
        MaxField = self.get_parameters()['field_setpoint']['maxval']
        MinField = self.get_parameters()['field_setpoint']['minval']
        MaxFieldSweep = self.get_parameters()['sweeprate_field']['maxval']
        MinFieldSweep = self.get_parameters()['sweeprate_field']['minval']
        # A to B conversion
        ABconversion = 115.733 / 14  # Ampere per Tesla
        self.add_parameter(
            'current_setpoint',
            type=types.FloatType,
            flags=Instrument.FLAG_GETSET | Instrument.FLAG_GET_AFTER_SET,
            minval=ABconversion * MinField,
            maxval=ABconversion * MaxField)
        self.add_parameter(
            'sweeprate_current',
            type=types.FloatType,
            flags=Instrument.FLAG_GETSET | Instrument.FLAG_GET_AFTER_SET,
            minval=ABconversion * MinFieldSweep,
            maxval=ABconversion * MaxFieldSweep)

        self.add_parameter(
            'remote_status',
            type=types.IntType,
            flags=Instrument.FLAG_GETSET | Instrument.FLAG_GET_AFTER_SET,
            format_map={
                0: "Local and locked",
                1: "Remote and locked",
                2: "Local and unlocked",
                3: "Remote and unlocked",
                4: "Auto-run-down",
                5: "Auto-run-down",
                6: "Auto-run-down",
                7: "Auto-run-down"})
        self.add_parameter('current', type=types.FloatType,
                           flags=Instrument.FLAG_GET)
        self.add_parameter('magnet_current', type=types.FloatType,
                           flags=Instrument.FLAG_GET)
        self.add_parameter('field', type=types.FloatType,
                           flags=Instrument.FLAG_GET)
        self.add_parameter('persistent_current', type=types.FloatType,
                           flags=Instrument.FLAG_GET)
        self.add_parameter('persistent_field', type=types.FloatType,
                           flags=Instrument.FLAG_GET)

        # Add functions
        self.add_function('get_all')
        self.get_all()

    def get_all(self):
        '''
        Reads all implemented parameters from the instrument,
        and updates the wrapper.

        Input:
            None

        Output:
            None
        '''
        logging.info(__name__ + ' : reading all settings from instrument')
        self.get_remote_status()
        self.get_current()
        self.get_magnet_current()
        self.get_current_setpoint()
        self.get_sweeprate_current()
        self.get_field()
        self.get_field_setpoint()
        self.get_sweeprate_field()
        self.get_persistent_current()
        self.get_persistent_field()
        self.get_activity()
        self.get_switch_heater()
        self.get_mode2()

    # Functions
    def _execute(self, message):
        '''
        Write a command to the device

        Input:
            message (str) : write command for the device

        Output:
            None
        '''
        logging.info(
            __name__ +
            ' : Send the following command to the device: %s' %
            message)
        self._visainstrument.write('@%s%s' % (self._number, message))
        sleep(20e-3)  # wait for the device to be able to respond
        result = self._visainstrument.read()
        if result.find('?') >= 0:
            print "Error: Command %s not recognized" % message
        else:
            return result

    def identify(self):
        '''
        Identify the device

        Input:
            None

        Output:
            None
        '''
        logging.info(__name__ + ' : Identify the device')
        return self._execute('V')

    def examine(self):
        '''
        Examine the status of the device

        Input:
            None

        Output:
            None
        '''
        logging.info(__name__ + ' : Examine status')

        print 'System Status: '
        print self.get_system_status()

        print 'Activity: '
        print self.get_activity()

        print 'Local/Remote status: '
        print self.get_remote_status()

        print 'Switch heater: '
        print self.get_switch_heater()

        print 'Mode: '
        print self.get_mode()

        print 'Polarity: '
        print self.get_polarity()

    def remote(self):
        '''
        Set control to remote and unlocked

        Input:
            None

        Output:
            None
        '''
        logging.info(__name__ + ' : Set control to remote and unlocked')
        self.set_remote_status(3)

    def local(self):
        '''
        Set control to local and unlcoked

        Input:
            None
        Output:
            None
        '''
        logging.info(__name__ + ' : Set control to local and unlocked')
        self.set_remote_status(2)

    def _do_get_remote_status(self):
        '''
        Get remote control status

        Input:
            None

        Output:
            result(str) :
            "Local & locked",
            "Remote & locked",
            "Local & unlocked",
            "Remote & unlocked",
            "Auto-run-down",
            "Auto-run-down",
            "Auto-run-down",
            "Auto-run-down"
        '''
        logging.info(__name__ + ' : Get remote control status')
        result = self._execute('X')
        return int(result[6])

    def _do_set_remote_status(self, mode):
        '''
        Set remote control status.

        Input:
            mode(int) :
            0 : "Local and locked",
            1 : "Remote and locked" (not available),
            2 : "Local and unlocked",
            3 : "Remote and unlocked",

        Output:
            None
        '''
        status = {
            0: "Local and locked",
            2: "Local and unlocked",
            3: "Remote and unlocked",
        }
        if status.__contains__(mode):
            logging.info(
                __name__ +
                ' : Setting remote control status to %s' %
                status.get(
                    mode,
                    "Unknown"))
            self._execute('C%s' % mode)
        else:
            print 'Invalid mode inserted: %s' % mode

    def get_system_status(self):
        '''
        Get the system status

        Input:
            None

        Output:
            result (str) :
            "Normal",
            "Quenched",
            "Over Heated",
            "Warming Up",
            "Fault"
        '''
        result = self._execute('X')
        logging.info(__name__ + ' : Getting system status')
        return int(result[1])

    def get_system_status2(self):
        '''
        Get the system status

        Input:
            None

        Output:
            result (str) :
            "Normal",
            "On positive voltage limit",
            "On negative voltage limit",
            "Outside negative current limit",
            "Outside positive current limit"
        '''
        result = self._execute('X')
        logging.info(__name__ + ' : Getting system status')
        return int(result[2])

    def _do_get_current(self):
        '''
        Demand output current of device

        Input:
            None

        Output:
            result (float) : output current in Amp
        '''
        logging.info(__name__ + ' : Read output current')
        result = self._execute('R0')
        return float(result.replace('R', ''))

    def get_voltage(self):
        '''
        Demand measured output voltage of device

        Input:
            None

        Output:
            result (float) : output voltage in Volt
        '''
        logging.info(__name__ + ' : Read output voltage')
        result = self._execute('R1')
        return float(result.replace('R', ''))

    def _do_get_magnet_current(self):
        '''
        Demand measured magnet current of device

        Input:
            None

        Output:
            result (float) : measured magnet current in Amp
        '''
        logging.info(__name__ + ' : Read measured magnet current')
        result = self._execute('R2')
        return float(result.replace('R', ''))

    def _do_get_current_setpoint(self):
        '''
        Return the set point (target current)

        Input:
            None

        Output:
            result (float) : Target current in Amp
        '''
        logging.info(__name__ + ' : Read set point (target current)')
        result = self._execute('R5')
        return float(result.replace('R', ''))

    def _do_set_current_setpoint(self, current):
        '''
        Set current setpoint (target current)
        Input:
            current (float) : target current in Amp

        Output:
            None
        '''
        logging.info(__name__ + ' : Setting target current to %s' % current)
        self.remote()
        self._execute('I%s' % current)
        self.local()
        self.get_field_setpoint()  # Update field setpoint

    def _do_get_sweeprate_current(self):
        '''
        Return sweep rate (current)

        Input:
            None

        Output:
            result (float) : sweep rate in Amp/min
        '''
        logging.info(__name__ + ' : Read sweep rate (current)')
        result = self._execute('R6')
        return float(result.replace('R', ''))

    def _do_set_sweeprate_current(self, sweeprate):
        '''
        Set sweep rate (current)

        Input:
            sweeprate(float) : Sweep rate in Amps/min.

        Output:
            None
        '''
        self.remote()
        logging.info(
            __name__ +
            ' : Set sweep rate (current) to %s Amps/min' %
            sweeprate)
        self._execute('S%s' % sweeprate)
        self.local()
        self.get_sweeprate_field()  # Update sweeprate_field

    def _do_get_field(self):
        '''
        Demand output field

        Input:
            None

        Output:
            result (float) : Output field in Tesla
        '''
        logging.info(__name__ + ' : Read output field')
        result = self._execute('R7')
        return float(result.replace('R', ''))

    def _do_get_field_setpoint(self):
        '''
        Return the set point (target field)

        Input:
            None

        Output:
            result (float) : Field set point in Tesla
        '''
        logging.info(__name__ + ' : Read field set point')
        result = self._execute('R8')
        return float(result.replace('R', ''))

    def _do_set_field_setpoint(self, field):
        '''
        Set the field set point (target field)
        Input:
            field (float) : target field in Tesla

        Output:
            None
        '''
        logging.info(__name__ + ' : Setting target field to %s' % field)
        self.remote()
        self._execute('J%s' % field)
        self.local()
        self.get_current_setpoint()  # Update current setpoint

    def _do_get_sweeprate_field(self):
        '''
        Return sweep rate (field)

        Input:
            None

        Output:
            result (float) : sweep rate in Tesla/min
        '''
        logging.info(__name__ + ' : Read sweep rate (field)')
        result = self._execute('R9')
        return float(result.replace('R', ''))

    def _do_set_sweeprate_field(self, sweeprate):
        '''
        Set sweep rate (field)

        Input:
            sweep rate in T/min
            sweeprate(float) : Sweep rate in Tesla/min.

        Output:
            None
        '''
        logging.info(
            __name__ +
            ' : Set sweep rate (field) to %s Tesla/min' %
            sweeprate)
        self.remote()
        self._execute('T%s' % sweeprate)
        self.local()
        self.get_sweeprate_current()  # Update sweeprate_current

    def get_voltage_limit(self):
        '''
        Return voltage limit

        Input:
            None

        Output:
            result (float) : voltage limit in Volt
        '''
        logging.info(__name__ + ' : Read voltage limit')
        result = self._execute('R15')
        result = float(result.replace('R', ''))
        self.set_parameter_bounds('voltage', -result, result)
        return result

    def _do_get_persistent_current(self):
        '''
        Return persistent magnet current

        Input:
            None

        Output:
            result (float) : persistent magnet current in Amp
        '''
        logging.info(__name__ + ' : Read persistent magnet current')
        result = self._execute('R16')
        return float(result.replace('R', ''))

    def get_trip_current(self):
        '''
        Return trip current

        Input:
            None

        Output:
            result (float) : trip current om Amp
        '''
        logging.info(__name__ + ' : Read trip current')
        result = self._execute('R17')
        return float(result.replace('R', ''))

    def _do_get_persistent_field(self):
        '''
        Return persistent magnet field

        Input:
            None

        Output:
            result (float) : persistent magnet field in Tesla
        '''
        logging.info(__name__ + ' : Read persistent magnet field')
        result = self._execute('R18')
        return float(result.replace('R', ''))

    def get_trip_field(self):
        '''
        Return trip field

        Input:
            None

        Output:
            result (float) : trip field in Tesla
        '''
        logging.info(__name__ + ' : Read trip field')
        result = self._execute('R19')
        return float(result.replace('R', ''))

    def get_heater_current(self):
        '''
        Return switch heater current

        Input:
            None

        Output:
            result (float) : switch heater current in milliAmp
        '''
        logging.info(__name__ + ' : Read switch heater current')
        result = self._execute('R20')
        return float(result.replace('R', ''))

    def get_current_limit_upper(self):
        '''
        Return safe current limit, most positive

        Input:
            None

        Output:
            result (float) : safe current limit, most positive in Amp
        '''
        logging.info(__name__ + ' : Read safe current limit, most positive')
        result = self._execute('R22')
        return float(result.replace('R', ''))

    def get_current_limit_lower(self):
        '''
        Return safe current limit, most negative

        Input:
            None

        Output:
            result (float) : safe current limit, most negative in Amp
        '''
        logging.info(__name__ + ' : Read safe current limit, most negative')
        result = self._execute('R21')
        return float(result.replace('R', ''))

    def get_lead_resistance(self):
        '''
        Return lead resistance

        Input:
            None

        Output:
            result (float) : lead resistance in milliOhm
        '''
        logging.info(__name__ + ' : Read lead resistance')
        result = self._execute('R23')
        return float(result.replace('R', ''))

    def get_magnet_inductance(self):
        '''
        Return magnet inductance

        Input:
            None

        Output:
            result (float) : magnet inductance in Henry
        '''
        logging.info(__name__ + ' : Read magnet inductance')
        result = self._execute('R24')
        return float(result.replace('R', ''))

    def _do_get_activity(self):
        '''
        Get the activity of the magnet. Possibilities: Hold, Set point, Zero or Clamp.
        Input:
            None

        Output:
            result(str) : "Hold", "Set point", "Zero" or "Clamp".
        '''
        result = self._execute('X')
        logging.info(__name__ + ' : Get activity of the magnet.')
        return int(result[4])

    def _do_set_activity(self, mode):
        '''
        Set the activity to Hold, To Set point or To Zero.

        Input:
            mode (int) :
            0 : "Hold",
            1 : "To set point",
            2 : "To zero"

            4 : "Clamped" (not included)

        Output:
            None
        '''
        status = {
            0: "Hold",
            1: "To set point",
            2: "To zero"
        }
        if status.__contains__(mode):
            logging.info(
                __name__ +
                ' : Setting magnet activity to %s' %
                status.get(
                    mode,
                    "Unknown"))
            self.remote()
            self._execute('A%s' % mode)
            self.local()
        else:
            print 'Invalid mode inserted.'

    def hold(self):
        '''
        Set the device activity to "Hold"
        Input:
            None
        Output:
            None
        '''
        self.set_activity(0)

    def to_setpoint(self):
        '''
        Set the device activity to "To set point". This initiates a sweep.
        Input:
            None
        Output:
            None
        '''
        self.set_activity(1)

    def to_zero(self):
        '''
        Set the device activity to "To zero". This sweeps te magnet back to zero.
        Input:
            None
        Output:
            None
        '''
        self.set_activity(2)

    def _do_get_switch_heater(self):
        '''
        Get the switch heater status.
        Input:
            None

        Output:
            result(str) : "Off magnet at zero", "On (switch open)", "Off magnet at field (switch closed)",
            "Heater fault (heater is on but current is low)" or "No switch fitted".
        '''
        logging.info(__name__ + ' : Get switch heater status')
        result = self._execute('X')
        return int(result[8])

    def _do_set_switch_heater(self, mode):
        '''
        Set the switch heater Off or On. Note: After issuing a command it is necessary to wait
        several seconds for the switch to respond.
        Input:
            mode (int) :
            0 : "Off",
            1 : "On, if PSU = Magnet",

            2 : "On, No checks" (not available)

        Output:
            None
        '''
        status = {
            0: "Off",
            1: "On, if PSU = Magnet"
        }
        if status.__contains__(mode):
            logging.info(
                __name__ +
                ' : Setting switch heater to %s' %
                status.get(
                    mode,
                    "Unknown"))
            self.remote()
            self._execute('H%s' % mode)
            print "Setting switch heater... (wait 40s)"
            self.local()
            sleep(40)
        else:
            print 'Invalid mode inserted.'
        sleep(0.1)
        self.get_switch_heater()

    def heater_on(self):
        '''
        Switch the heater on, with PSU = Magnet current check
        Input:
            None
        Output:
            None
        '''

        current_in_magnet = self.get_persistent_current()
        current_in_leads = self.get_current()
        if self.get_switch_heater() == 1:
            print 'Heater is already on!'
        else:
            if self.get_mode2() == 0:

                if current_in_leads == current_in_magnet:

                    self.set_switch_heater(1)
                else:
                    print 'Current in the leads is not matching persistent current!'
            else:
                print 'Magnet supply not at rest, cannot switch on heater!'

        self.get_switch_heater()

    def set_persistent(self):
        '''
        Puts magnet into persistent mode
        '''

        if self.get_mode() == 0:
            self.heater_off()
            sleep(20)
            self.to_zero()
            self.get_all()
        else:
            print 'Magnet is not at rest, cannot put it in persistent mode'
        self.get_all()

    def leave_persistent_mode(self):
        '''
        Read out persistent current, match the current in the leads to that current
        and switch on heater
        Input:
            None
        Output:
            None
        '''
        if self.get_switch_heater() == 2:
            current_in_magnet = self.get_persistent_current()
            current_in_leads = self.get_current()
            self.hold()
            self.set_current_setpoint(current_in_magnet)
            self.to_setpoint()

            while current_in_leads != current_in_magnet:
                current_in_leads = self.get_current()
            self.heater_on()
            self.hold()

        elif self.get_switch_heater() == 1:
            print 'Heater is already on, so the magnet was not in persistent mode'
        elif self.get_switch_heater() == 0:
            print 'Heater is off, but magnet is not in persistent mode. Please, check magnet locally!'

        self.get_all()
        self.get_changed()

    def run_to_field(self, field_value):
        if self.get_switch_heater() == 1:
            self.hold()
            self.set_field_setpoint(field_value)
            self.to_setpoint()
        else:
            print 'Switch heater is off, cannot change the field.'
        self.get_all()

    def run_to_field_wait(self, field_value):
        '''
        Please comment this...
        Go to field_value and wait until it's done sweeping.
        '''
        if self.get_switch_heater() == 1:
            self.hold()
            self.set_field_setpoint(field_value)
            self.to_setpoint()
            self.local()
            magnet_mode = self.get_mode2()
            while magnet_mode != 0:
                print 'Magnet still sweeping, current field %s' % self.get_field()
                magnet_mode = self.get_mode2()
                sleep(0.5)
        else:
            print 'Switch heater is off, cannot change the field.'
        self.get_all()
        self.local()

    def heater_off(self):
        '''
        Switch the heater off
        Input:
            None
        Output:
            None
        '''
        if self.get_switch_heater() == 0 | 2:
            print 'Heater is already off!'
        else:

            if self.get_mode2() == 0:

                self.set_switch_heater(0)
            else:
                print 'Magnet is not at rest, cannot switch of the heater!'

    def get_mode(self):
        '''
        Get the Mode of the device
        Input:
            None

        Output:
            "Amps, Magnet sweep: fast",
            "Tesla, Magnet sweep: fast",
            "Amps, Magnet sweep: slow",
            "Tesla, Magnet sweep: slow"
        '''
        logging.info(__name__ + ' : Get device mode')
        result = self._execute('X')
        return int(result[10])

    def _do_get_mode2(self):
        '''
        Get the Mode of the device
        Input:
            None

        Output:
            "At rest",
            "Sweeping",
            "Sweep limiting",
            "Sweeping & sweep limiting"
        '''
        logging.info(__name__ + ' : Get device mode')
        result = self._execute('X')
        return int(result[11])

    def set_mode(self, mode):
        '''
        Input:
            mode(int):
            0 : "Amps, Magnet sweep: fast",
            1 : "Tesla, Magnet sweep: fast",
            4 : "Amps, Magnet sweep: slow",
            5 : "Tesla, Magnet sweep: slow"
            8 : "Amps, (Magnet sweep: unaffected)",
            9 : "Tesla, (Magnet sweep: unaffected)"

        Output:
            None
        '''
        status = {
            0: "Amps, Magnet sweep: fast",
            1: "Tesla, Magnet sweep: fast",
            4: "Amps, Magnet sweep: slow",
            5: "Tesla, Magnet sweep: slow",
            8: "Amps, (Magnet sweep: unaffected)",
            9: "Tesla, (Magnet sweep: unaffected)"
        }
        if status.__contains__(mode):
            logging.info(
                __name__ +
                ' : Setting device mode to %s' %
                status.get(
                    mode,
                    "Unknown"))
            self.remote()
            self._execute('M%s' % mode)
            self.local()
        else:
            print 'Invalid mode inserted.'

    def get_polarity(self):
        '''
        Get the polarity of the output current
        Input:
            None

        Output:
            result (str) :
            "Desired: Positive, Magnet: Positive, Commanded: Positive",
            "Desired: Positive, Magnet: Positive, Commanded: Negative",
            "Desired: Positive, Magnet: Negative, Commanded: Positive",
            "Desired: Positive, Magnet: Negative, Commanded: Negative",
            "Desired: Negative, Magnet: Positive, Commanded: Positive",
            "Desired: Negative, Magnet: Positive, Commanded: Negative",
            "Desired: Negative, Magnet: Negative, Commanded: Positive",
            "Desired: Negative, Magnet: Negative, Commanded: Negative"
        '''
        status1 = {
            0: "Desired: Positive, Magnet: Positive, Commanded: Positive",
            1: "Desired: Positive, Magnet: Positive, Commanded: Negative",
            2: "Desired: Positive, Magnet: Negative, Commanded: Positive",
            3: "Desired: Positive, Magnet: Negative, Commanded: Negative",
            4: "Desired: Negative, Magnet: Positive, Commanded: Positive",
            5: "Desired: Negative, Magnet: Positive, Commanded: Negative",
            6: "Desired: Negative, Magnet: Negative, Commanded: Positive",
            7: "Desired: Negative, Magnet: Negative, Commanded: Negative"
        }
        status2 = {
            1: "Negative contactor closed",
            2: "Positive contactor closed",
            3: "Both contactors open",
            4: "Both contactors closed"
        }
        logging.info(__name__ + ' : Get device polarity')
        result = self._execute('X')
        return status1.get(int(result[13]), "Unknown") + \
            ", " + status2.get(int(result[14]), "Unknown")

    def get_changed(self):
        print "Current: "
        print self.get_current()
        print "Field: "
        print self.get_field()
        print "Magnet current: "
        print self.get_magnet_current()
        print "Heater current: "
        print self.get_heater_current()
        print "Mode: "
        print self.get_mode()
