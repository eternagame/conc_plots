import os, random, string, glob, itertools
import numpy as np
import subprocess as sp

NUPACKDIR = '%s/bin' % os.environ['NUPACKHOME']
TMPDIR = '/tmp'

def _filename(n=6):
    """generate random filename

    Args:
        n (int): number of characters
    """
    rand = ''.join([random.choice(string.ascii_lowercase) for _ in range(n)])
    return '%s/%s' % (TMPDIR, rand)

def _write(lines, fname=None):
    """write lines to file

    Args:
        lines (list): line(s) to write to file
        fname (str): filename to write to
    """
    if fname is None:
        fname = '%s.in' % _filename()
    with open(fname, 'w') as f:
        for line in lines:
            f.write('%s\n' % line)
    return fname

def _parse_ppairs_chunk(f):
    """parse one set of bpps from ppairs

    Args:
        f (file): open ppairs file object

    Returns:
        str, np.ndarray: complex identifier and array containing base pair
                         probabilities
    """
    key = f.readline().strip()[2:]
    n = int(f.readline())
    bpp = np.zeros((n, n))
    line = f.readline()
    while not line.startswith('% %'):
        s = line.split()
        if int(s[1]) != n + 1:
            bpp[int(s[0])-1, int(s[1])-1] = float(s[2])
        line = f.readline()
    return key, bpp

def _parse_ppairs(fname):
    """parse base pair probability file

    Args:
        fname (str): path of ppairs output file

    Returns:
        dict: key is complex identifier, value is np.ndarray containing
              base pair probabilities
    """
    bpps = {}
    with open(fname) as f:
        line = f.readline()
        while line:
            if line.startswith('% %'):
                key, value = _parse_ppairs_chunk(f)
                bpps[key] = value
            line = f.readline()
    return bpps

def _get_orders(n):
    """get all possible orders of n strands"""
    all_ = []
    # loop over number of strands
    for i in range(1,n+1):
        # loop over each possible combination
        for strands in itertools.combinations(range(1,n+1), i):
            if len(strands) > 2:
                # loop over all orders, fixing first strand to avoid circular
                # permutations
                for order in itertools.permutations(strands[1:]):
                    all_.append([strands[0]] + list(order))
            else:
                all_.append(list(strands))
    return all_

def complexes(seqs, order=None, molecule='dna', T=37, keep=False, bpp=False):
    """get free energy of multi-strand complexes

    Args:
        seqs (str or list): input filename or list of nucleic acid sequences
        order (list): list of strand orderings
        molecule (str): "dna" or "rna"
        T (float): temperature
        keep (bool): whether or not to keep output files
        bpp (bool): whether or not to compute base pair probabilities

    Returns:
        dict: each key, value pair represents one species; the key is a list
              representing the counts for each strand, value is the free energy
              and bpp matrix if bpp is true
    """
    if isinstance(seqs, str):
        fname = seqs
        nseqs = int(sp.check_output(['wc', '-l', '%s.in' % seqs]).split()[0]) - 2
    else:
        fname = _write([str(len(seqs))] + seqs + [1])[:-3]
        nseqs = len(seqs)
    command = ['%s/complexes' % NUPACKDIR, '-material', molecule, '-T', str(T),
               '-ordered']
    if bpp:
        command.append('-pairs')
    if order is None:
        order = _get_orders(nseqs)
    _write([' '.join([str(s) for s in strands]) for strands in order],
          '%s.list' % fname)
    p = sp.Popen(command + [fname], stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode:
        raise Exception('complexes failed: on %s\n%s' % (seqs, stderr))
    if bpp:
        bpps = _parse_ppairs('%s.ocx-ppairs' % fname)
    strands = []
    with open('%s.ocx-key' % fname) as f:
        line = f.readline()
        while line:
            if not line.startswith('%'):
                strands.append(tuple([int(x) for x in line.split()[2:]]))
            line = f.readline()
    with open('%s.ocx' % fname) as f:
        line = f.readline()
        while not line.strip() or line.startswith('%'):
            line = f.readline()
        result = {}
        while line:
            spl = line.split()
            try:
                energy = float(spl[-1])
            except:
                line = f.readline()
                continue
            if bpp:
                result[strands.pop(0)] = [
                    energy,
                    bpps['complex%s-order%s' % (spl[0], spl[1])]]
            else:
                result[strands.pop(0)] = energy
            line = f.readline()
    if not keep:
        for f in glob.glob('%s*' % fname):
            os.remove(f)
    return result
