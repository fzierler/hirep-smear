#!/usr/bin/env python3
'''
This script converts the occurrences of '.im' and '.re' into 
usages of _complex_im and _complex_re.

WARNING: This is not what is needed for lvalues.
'''
import re
from sys import argv
text = open(argv[1],'r').read()

# patterns of the kind 
# (u).$cname\[$idx\].re 

text1 = re.subn('(?P<arg>\(.\)\.[_$a-zA-Z0-9\\\[\]]+)\.re','_complex_re(\g<arg>)', text )[0]

# (u).$cname\[$idx\].im
text2 = re.subn('(?P<arg>\(.\)\.[_$a-zA-Z0-9\\\[\]]+)\.im','_complex_im(\g<arg>)', text1)[0]

# _tmp[6].re
text3 = re.subn('(?P<arg>_[a-zA-Z_]+\[[0-9]+\])\.re','_complex_re(\g<arg>)', text2)[0]

# _tmp[6].im
textfinal = re.subn('(?P<arg>_[a-zA-Z_]+\[[0-9]+\])\.im','_complex_im(\g<arg>)', text3)[0]


print(textfinal)
