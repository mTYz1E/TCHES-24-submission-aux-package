import chipwhisperer as cw
import math
import datetime
import time
import os
import sys
import subprocess

# run bash command and return output
def run_cmd(cmd):
    return subprocess.check_output(cmd, shell=True).decode('utf-8')

NSHARES = 4
DECIMATE=1
TRACE_NUM = 500
TRACE_FILE = 'masked-cmp-{}.bin'.format(datetime.datetime.now().strftime('%Y-%m-%d-%H:%M:%S'))

try:
    if not scope.connectStatus:
        scope.con()
except NameError:
    scope = cw.scope()

PLATFORM='CW308_STM32F4'
CRYPTO_TARGET = 'NONE'
SS_VER = 'SS_VER_2_1'

try:
    if SS_VER == "SS_VER_2_1":
        target_type = cw.targets.SimpleSerial2
    elif SS_VER == "SS_VER_2_0":
        raise OSError("SS_VER_2_0 is deprecated. Use SS_VER_2_1")
    else:
        target_type = cw.targets.SimpleSerial
except:
    SS_VER="SS_VER_1_1"
    target_type = cw.targets.SimpleSerial

try:
    target = cw.target(scope, target_type)
except:
    print("INFO: Caught exception on reconnecting to target - attempting to reconnect to scope first.")
    print("INFO: This is a work-around when USB has died without Python knowing. Ignore errors above this line.")
    scope = cw.scope()
    target = cw.target(scope, target_type)


print("INFO: Found ChipWhisperer")

if "STM" in PLATFORM or PLATFORM == "CWLITEARM" or PLATFORM == "CWNANO":
    prog = cw.programmers.STM32FProgrammer
elif PLATFORM == "CW303" or PLATFORM == "CWLITEXMEGA":
    prog = cw.programmers.XMEGAProgrammer
elif "neorv32" in PLATFORM.lower():
    prog = cw.programmers.NEORV32Programmer
elif PLATFORM == "CW308_SAM4S":
    prog = cw.programmers.SAM4SProgrammer
else:
    prog = None

TRACE_SAMPLES_NUM = 24400

time.sleep(0.05)
scope.default_setup()
# overwrite some defaults
# use down sampling to record the whole trace
scope.adc.decimate = DECIMATE
scope.adc.timeout = 2
scope.adc.samples = TRACE_SAMPLES_NUM

print(scope.adc)

def reset_target(scope):
    if PLATFORM == "CW303" or PLATFORM == "CWLITEXMEGA":
        scope.io.pdic = 'low'
        time.sleep(0.1)
        scope.io.pdic = 'high_z' #XMEGA doesn't like pdic driven high
        time.sleep(0.1) #xmega needs more startup time
    elif "neorv32" in PLATFORM.lower():
        raise IOError("Default iCE40 neorv32 build does not have external reset - reprogram device to reset")
    elif PLATFORM == "CW308_SAM4S":
        scope.io.nrst = 'low'
        time.sleep(0.25)
        scope.io.nrst = 'high_z'
        time.sleep(0.25)
    else:
        scope.io.nrst = 'low'
        time.sleep(0.05)
        scope.io.nrst = 'high_z'
        time.sleep(0.05)

# flash the firmware
fw_path = './kyber-masked-cmp-{}.hex'.format(PLATFORM)
cw.program_target(scope, prog, fw_path)

time.sleep(5)

target.flush()

def simpleserial_get_bin(target):
    reply = target.read_cmd('r')
    data_len = reply[2]
    data = reply[3:3 + data_len]
    return data

match NSHARES:
    case 2:
        WAIT_SEC = 0.2
    case 3:
        WAIT_SEC = 0.7
    case 4:
        WAIT_SEC = 1.3
    case _:
        print('untested masking order: {}'.format(NSHARES))
        sys.exit(0)

print('number of traces: {}'.format(TRACE_NUM))
print('masking order: {}'.format(NSHARES))
print('interval between measurements: {}'.format(WAIT_SEC))
print('decimate: {}'.format(scope.adc.decimate))

print('start collecting traces')
with open(TRACE_FILE, 'wb') as output_f:
    output_f.write(TRACE_NUM.to_bytes(4))
    output_f.write(TRACE_SAMPLES_NUM.to_bytes(4))

    for i in range(TRACE_NUM):
        target.send_cmd('n', 0, bytearray())
        target.simpleserial_wait_ack()

        scope.arm()
        target.send_cmd('c', 0, bytearray())

        time.sleep(WAIT_SEC)
        target.simpleserial_wait_ack()

        ret = scope.capture()
        if ret:
            print('target time out')
            continue

        # NOTE: A trace is an array of type nympy.float64
        tr_no_dec = scope.get_last_trace()
        if len(tr_no_dec) != TRACE_SAMPLES_NUM:
            print('samples in the trace less than expected: {}'.format(len(tr_no_dec)))
            continue

        output_f.write(tr_no_dec)

        target.send_cmd('f', 0, bytearray())
        target.simpleserial_wait_ack()

        scope.arm()
        target.send_cmd('c', 0, bytearray())

        time.sleep(WAIT_SEC)
        target.simpleserial_wait_ack()

        ret = scope.capture()
        if ret:
            print('target time out')
            continue

        tr_dec = scope.get_last_trace()
        if len(tr_dec) != TRACE_SAMPLES_NUM:
            print('samples in the trace less than expected: {}'.format(len(tr_dec)))
            continue

        output_f.write(tr_dec)
