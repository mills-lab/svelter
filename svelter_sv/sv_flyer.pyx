#!/usr/bin/env python

from __future__ import print_function
import re
import os
import string
import pandas as pd
import itertools
ttab = string.maketrans(string.ascii_uppercase[:], string.ascii_lowercase[::-1])
inv_dict = {
    "a^" : "Z", "b^" : "Y", "c^" : "X", "d^" : "W", "e^" : "V",
    "f^" : "U", "g^" : "T", "h^" : "S", "i^" : "R", "j^" : "Q",
    "k^" : "P", "l^" : "O", "m^" : "N", "n^" : "M", "o^" : "L",
    "p^" : "K", "q^" : "J", "r^" : "I", "s^" : "H", "t^" : "G",
    "u^" : "F", "v^" : "E", "w^" : "D", "x^" : "C", "y^" : "B",
    "z^" : "A"}

# sv_array schema:      [dlt,    inv,   tan,    dis,   inv_tan, inv_dis, transl, transp]
templates = pd.DataFrame({
        'del_inv' :     [True,  True,  False, False, False,  False,  False, False],
        'tan_dup_inv' : [False, True,  False, False, True,   False,  False, False],
        'dis_dup_inv' : [False, True, False, False, False,  True,   False, False],
        'del_tan_dup' : [True,  False, True,  False, False,  False,  False, False],
        'del_dis_dup' : [True,  False, False, True,  False,  False,  False, False],
        'transl' :      [False, False, False, False, False,  False,  True,  False],
        'transl_inv' :  [False, True, False, False, False,  False,  True,  False],
        'transp' :      [False, False, False, False, False,  False,  False, True],
        'dis_dup' :     [False, False, False, True,  False,  False,  False, False],
        'tan_dup' :     [False, False, True,  False, False,  False,  False, False],
        'com_dup' :     [False, False, True,  True,  False,  False,  False, False],
        'inv' :         [False, True,  False, False, False,  False,  False, False],
        'del' :         [True,  False, False, False, False,  False,  False, False],
        'undeter' :     [False, False, False, False, False, False, False, False],
        'nc' :          [None,  None,  None,  None,   None,  None,  None,  None]
        # 'complex':    combination of above arrays
        })


def inverse_shadow(STR):
    """converts string inversions to capital letter analogs
    
    Args:
        STR (str): the haplotype with inversion to be converted
    
    Returns:
        inv (str): converted string
    """
    global inv_dict
    matches = re.finditer(r"\^", STR)
    idx = [(match.start()-1, match.end()) for match in matches]
    subset = set((STR[span[0]:span[1]], inv_dict.get(STR[span[0]:span[1]])) for span in idx[::-1])
    inv = STR[:]
    for sub in subset:
        inv = inv.replace(sub[0], sub[1])
    return inv
    
    
def knee_pads(STR):
    """Generator that creates a series of sliding windows
    
    Args:
        STR (str): the sequence to be placed in sliding window
    
    Yeilds:
        an iteratively larger sliding window across the string until it is one less than the input
    """
    for N in reversed(range(len(STR))):
        if 1 < N < len(STR):
            i = iter(STR)
            win = tuple(itertools.islice(i, N))
            out = "".join(win)
            if len(win) == N:
                yield out
            for e in i:
                win = win[1:] + (e,)
                out = "".join(win)
                yield out
        else:
            pass
    
    
def clown_cartridge(ALT):
    """Checks for strange duplciation events
        e.g.
            ref = 'abc'
            alt = 'abcbc'
        This is equivalent to:
            ref = 'ad'
            alt = 'add'
    Args:
        ALT (str): alternate haplotype
    
    Returns:
        Boolean/string of strange duplication event [True = Tandem, False = Disperse, None = None, Complex = "Complex"])
    """
    slider = knee_pads(ALT)
    tanc = set()
    for i in slider:
        win_len = len(i)
        if ALT.count(i) > 1:
            tmp = [(e.start(), e.end()) for e in list(re.finditer(i, ALT))]
            if all(tmp[j+1][0] == tmp[j][0] + win_len for j in range(len(tmp)-1)):
                tanc.update((True,))
            elif any(tmp[j+1][0] == tmp[j][0] + win_len for j in range(len(tmp)-1)):
                tanc.update((True, False))
            else:
                tanc.update((False,))
        else:
            pass
        if (True or False) in tanc:
            break
    if all(i == True for i in tanc):
        return True
    elif all(i == False for i in tanc):
        return False
    elif len(tanc) == 0:
        return None
    else:
        return "Complex"
    
    
def perfect_balance(REF, ALT):
    """Compares the reference to alternate HAPLOTYPE
    
    Args:
        REF (str): reference haplotype
        ALT (str); alternate haplotype
    
    Returns:
        a boolean array representing ontological description
    """
    global ttab
    inverted = ALT[:].translate(ttab)
    dlt, inv, tan, dis, inv_tan, inv_dis, transl, transp = False, False, False, False, False, False, False, False
    
    if ALT == REF:
        same = True
    else:
        same = False
    
    if not same:
        # check to see if deletion occured
        dlt = (True if any([len(set(REF).difference(set(inverted))) >= 1,
                not any(letter.isupper() for letter in ALT) and len(
                set(REF).difference(set(ALT))) >=1]) else False)
        
        # check for inversions
        inv = (True if any(letter.isupper() for letter in ALT) else False)
        
        # check for transposable elements
        transp = (True if sum(1 for letter in set(ALT) if letter not 
                in REF and not letter.isupper()) >= 1 else False)
        
        # check for translocation
        transl = (True if all([any(ALT[i] > ALT[i+1] for i in range(len(ALT)-1)),
                len(REF) == len(ALT), set(REF) == set(ALT)]) else False)
        
        # check for strangely defined tandem or disperse duplications
        bad_dup = clown_cartridge(ALT)
        if bad_dup is True:
            tan = True
            dis = False
        elif bad_dup is False:
            tan = False
            dis = True
        else:
            # check for simple tandem or disperse duplications
            tan, dis = False, False
            while not (tan and dis):
                for i in range(len(ALT)-1):
                    tan_check = False
                    if ALT[i] == ALT[i+1]:
                        tan_check = True
                        tan = True
                    if ALT[i] in ALT[i+2:] and not tan_check:
                        dis = True
                break
        
        # check for inverted tandem, disperse duplication, or inverted translocation
        if inv:
            inv_bad_dup = clown_cartridge(inverted)
            transl = (True if all([any(inverted[i] > inverted[i+1] for i in range(len(inverted)-1)),
                len(REF) == len(inverted), set(REF) == set(inverted)]) else False)
            inv_tan, inv_dis = False, False
            if inv_bad_dup is True:
                inv_tan = True
                inv_dis = False
            elif inv_bad_dup is False:
                inv_tan = False
                inv_dis = True
            else:
                while not (inv_tan and inv_dis):
                    for i in range(len(inverted)-1):
                        tan_check = False
                        if inverted[i] == inverted[i+1]:
                            tan_check = True
                            inv_tan = True
                        if inverted[i] in inverted[i+2:] and not tan_check:
                            inv_dis = True
                    break
        return pd.Series([dlt, inv, tan, dis, inv_tan, inv_dis, transl, transp])
    else:
        return pd.Series([None, None, None, None, None, None, None, None])


def angry_magic(SV_ARRAY):
    global templates
    for name in templates:
        if (SV_ARRAY == templates[name]).all() == True:
            return name


def luck_in_the_chamber(ROW, IDX = None):
    """Pre-processes SVelter row
    
    Args:
        ROW (list): list containing the separated row from tab delimited file
        IDX (int): the index of the alternate allele prediction in list (if given)
        
    Returns:
        string value of matching classification
    
    """
    global templates
    one0 = None
    output = []
    if IDX is not None:
        ref = ROW[IDX-1].split('/')[0]
        alt = ROW[IDX].split('/')
    else:
        ref = ROW[4].split('/')[0]
        alt = ROW[5].split('/')
    
    if len(alt[0]) >= 1:
        if alt[0] == ref:
            output.append('nc')
        else:
            haplo = inverse_shadow(alt[0])
            booly = perfect_balance(ref, haplo)
            output.append(angry_magic(booly) if angry_magic(booly) is not None else 'complex')
    else:
        output.append("del")
    if len(alt[1]) >= 1:
        if alt[1] == ref:
            output.append('nc')
        else:
            haplo = inverse_shadow(alt[1])
            booly = perfect_balance(ref, haplo)
            output.append(angry_magic(booly) if angry_magic(booly) is not None else 'complex')
    else:
        output.append("del")
    return output


def classify(SVelter):
    """function for importing into SVelter
    
    Args:
        SVelter (str | list): the regular SVelter row output
    
    Returns:
        an appended row containing both alternate allele classifications as new columns
    """
    if type(SVelter) is str:
        row = SVelter.split('\t')
    elif type(SVelter) is list:
        row = SVelter[:]
    if 'start' in row[1]:
        row.extend(['alt1', 'alt2'])
        return row
    else:
        ROW = [item for item in row]
        flyer = luck_in_the_chamber(row)
        ROW.extend(flyer)
        return ROW

