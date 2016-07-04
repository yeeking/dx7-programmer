# DX7 Programmer: automatic DX7 programmer using the dexed DX7 emulator. 

This will render the presets from a synth and then try to find settings that make sounds like those presets:
```
cd python
python ga_find_test_set.py
```

You might need a couple of python packages:

```
pip install scipy
pip install scikits.talkbox
```

## External projects I am using:

* DX7 emulator: https://github.com/asb2m10/dexed
* Command line VST runner: https://github.com/teragonaudio/MrsWatson
* Other free VST plugins for testing: http://mda.smartelectronix.com
