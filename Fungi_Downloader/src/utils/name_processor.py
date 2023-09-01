import re


def process_two_part_name(name):
    name = name.lower().replace(' ', '_')
    name = name.replace('[', '')
    name = name.replace(']', '')
    name = re.sub(r'[^\w\s]+', '_', name)  # Replace non unicode chars with _
    name = re.sub(r'(_)\1+', '_', name)  # Replace consecutive __ with one _

    name = name.split('_')[:2]
    name = '_'.join(name)
    
    return name


def process_name(name):
    name = name.lower().replace(' ', '_')
    name = name.replace('[', '')
    name = name.replace(']', '')
    name = re.sub(r'[^\w\s]+', '_', name)  # Replace non unicode chars with _
    name = re.sub(r'(_)\1+', '_', name)  # Replace consecutive _
    
    return name