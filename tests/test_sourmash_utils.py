from sourmash_utils import *
from sourmash_utils import create_minhash_from_args, load_index_and_select
    
import os


def test_basic():
    pass


def test_minhash_create_dna():
    # DNA moltype works
    mh = FracMinHash(ksize=31, scaled=1000, moltype='DNA')
    assert mh.moltype == 'DNA'
    assert mh.scaled == 1000
    assert mh.ksize == 31


def test_minhash_create_dna_default():
    # DNA minhashes are default
    mh = FracMinHash(ksize=31, scaled=1000)
    assert mh.moltype == 'DNA'
    assert mh.scaled == 1000
    assert mh.ksize == 31


def test_minhash_create_protein():
    # test moltype protein
    mh = FracMinHash(ksize=7, scaled=100, moltype='protein')
    assert mh.moltype == 'protein'
    assert mh.scaled == 100
    assert mh.ksize == 7


def test_minhash_create_dayhoff():
    # test moltype dayhoff
    mh = FracMinHash(ksize=7, scaled=100, moltype='dayhoff')
    assert mh.moltype == 'dayhoff'
    assert mh.scaled == 100
    assert mh.ksize == 7


def test_minhash_create_hp():
    # test moltype hp
    mh = FracMinHash(ksize=7, scaled=100, moltype='hp')
    assert mh.moltype == 'hp'
    assert mh.scaled == 100
    assert mh.ksize == 7


def test_minhash_create_skipm1n3():
    # test moltype skipm1n3
    mh = FracMinHash(ksize=21, scaled=1000, moltype='skipm1n3')
    assert mh.moltype == 'skipm1n3'
    assert mh.scaled == 1000
    assert mh.ksize == 21


def test_minhash_create_skipm2n3():
    # test moltype skipm1n3
    mh = FracMinHash(ksize=21, scaled=1000, moltype='skipm2n3')
    assert mh.moltype == 'skipm2n3'
    assert mh.scaled == 1000
    assert mh.ksize == 21


def test_load_index_and_select():
    # test loading index and selecting without picklist
    select_mh = create_minhash_from_args(ksize=31,scaled=1000,moltype='DNA', dna=True, abund=True)

    here = os.path.dirname(__file__)
    fpath = os.path.join(here, 'podar-ref.zip')
    idx = load_index_and_select(filename=fpath, minhash_obj=select_mh)
    assert idx

def test_load_index_and_select_picklist():
    # test loading index and selecting with picklist
    select_mh = create_minhash_from_args(ksize=31,scaled=1000,moltype='DNA',abund=True)

    here = os.path.dirname(__file__)
    fpath = os.path.join(here, 'podar-ref.zip')
    ppath = os.path.join(here, 'podar-ref.picklist.csv')
    idx = load_index_and_select(filename=fpath, minhash_obj=select_mh, pickfile=ppath, coltype='ident', colname='accession')
    assert idx
    assert len(idx) == 9
