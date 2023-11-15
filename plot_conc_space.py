import argparse
import pandas as pd
import numpy as np
import nupack
import multiprocessing as mp
import itertools as it
import matplotlib as mpl
import matplotlib.pyplot as plt

rt = 0.616


def parse_states(filename):
    """
    parse information about inputs, reporter, and complexes to use in
    simulation

    Args:
        filename (str): name of file to read info from

    Returns:
        inputs (list): list of input sequences
        reporter (str): reporter sequence
        orders (dict): complexes to consider
    """
    inputs = {}
    orders = {}
    with open(filename) as f:
        reporter = f.readline().strip()
        inputs = {chr(i+65): s for i,s in enumerate(f.readline().split())}
        line = f.readline()
        while line:
            strands = [int(s) for s in line.split()]
            orders[''.join([chr(int(s)+63) for s in strands])] = strands
            line = f.readline()
    inputs = [inputs[k] for k in sorted(inputs)]
    return inputs, reporter, orders


def get_pfs(sequence, inputs, reporter, complexes):
    """
    get partition function information for one sequence

    Args:
        sequence (str): RNA sequence
        inputs (list): list of input sequences
        reporter (str): reporter sequence
        orders (dict): complexes to consider

    Returns:
        pd.Series: contains fields with free enrgies of each state
    """
    seqs = [sequence] + inputs + [reporter]
    rcomplexes = {k: [[1] + list(x)
                  for x in it.permutations(v + [len(seqs)], len(v) + 1)]
                  for k,v in complexes.items()}
    complexes = {k: [[1] + list(x) for x in it.permutations(v, len(v))]
                 for k,v in complexes.items()}
    allorders = [order for orders in rcomplexes.values() for order in orders] + \
                [order for orders in complexes.values() for order in orders] + \
                [(1,), (1, len(seqs))]
    result = nupack.complexes(seqs, molecule='rna', order=allorders)
    s = pd.Series()
    for name in complexes:
        orders = complexes[name]
        rorders = rcomplexes[name]
        dGs = np.full(len(orders), np.inf)
        rdGs = np.full(len(rcomplexes[name]), np.inf)
        for i, order in enumerate(orders):
            dGs[i] = result[tuple(order)]
        for i, order in enumerate(rorders):
            rdGs[i] = result[tuple(order)]
        s['dG_%s' % name] = -rt*np.log(np.sum(np.exp(-dGs/rt)))
        s['dG_%s_reporter' % name] = -rt*np.log(np.sum(np.exp(-rdGs/rt)))
    s['dG_none'] = result[(1,)]
    s['dG_none_reporter'] = result[(1, len(seqs))]
    return s


def get_pfs_multiple(seqs, inputs=None, reporter=None, complexes=None,
                     inputfile=None):
    """
    get partition function information for multiple sequences

    Args:
        sequence (list): list of RNA sequences
        inputs (list): list of input sequences
        reporter (str): reporter sequence
        orders (dict): complexes to consider

    Returns:
        pd.Series: contains fields with free enrgies of each state
    """
    if inputfile is None:
        if inputs is None or reporter is None or complexes is None:
            raise ValueError('if inputfile not provided, inputs, reporter, '
                             'and complexes must be provided')
    else:
        inputs, reporter, complexes = parse_states(inputfile)

    p = mp.Pool()
    results = p.map(get_pfs, it.izip(
        seqs, it.repeat(inputs), it.repeat(reporter),
        it.repeat(complexes)))
    return pd.concat(results, axis=1).T


def get_conc_space(s, nsteps=25, reporter=1e-9):
    """
    get simulation over 3d concentration space

    Args:
        s (pd.Series): contains fields with free energies of each state
        nsteps (int): number of increments along each axis
        reporter (float): concentration of reporter

    Returns:
        np.ndarray: tensor of proportion bound to reporter
    """
    states = [col for col in s.index if col.startswith('dG') and
              '%s_reporter' % col in s.index]
    concrange = np.logspace(-15, -3, nsteps)
    num = np.zeros((nsteps, nsteps, nsteps))
    denom = np.zeros((nsteps, nsteps, nsteps))
    conc = {'A': np.tile(concrange[:,np.newaxis,np.newaxis], (1,nsteps,nsteps)),
            'B': np.tile(concrange[np.newaxis,:,np.newaxis], (nsteps,1,nsteps)),
            'C': np.tile(concrange[np.newaxis,np.newaxis,:], (nsteps,nsteps,1))}
    for state in states:
        multiplier = np.ones((nsteps, nsteps, nsteps))
        if state != 'dG_none':
            for c in state[3:]:
                multiplier = multiplier * conc[c]
        num += np.exp(-s['%s_reporter' % state]/rt) * multiplier * reporter
        denom += np.exp(-s[state]/rt) * multiplier + \
            np.exp(-s['%s_reporter' % state]/rt) * multiplier * reporter
    return num/denom


def get_2d_space(s, nsteps=25, reporter=1e-9):
    """
    get simulation over 2d space - AB vs CC

    Args:
        s (pd.Series): contains fields with free energies of each state
        nsteps (int): number of increments along each axis
        reporter (float): concentration of reporter

    Returns:
        np.ndarray: tensor of proportion bound to reporter
    """
    states = ['dG_none', 'dG_AB', 'dG_CC', 'dG_ABCC']
    concrange = np.logspace(-30, -6, nsteps)
    num = np.zeros((nsteps, nsteps))
    denom = np.zeros((nsteps, nsteps))
    conc = {'AB': np.tile(concrange[:,np.newaxis], (1,nsteps)),
            'CC': np.tile(concrange[np.newaxis,:], (nsteps,1))}
    for state in states:
        multiplier = np.ones((nsteps, nsteps))
        if 'AB' in state:
            multiplier = multiplier * conc['AB']
        if 'CC' in state:
            multiplier = multiplier * conc['CC']
        num += np.exp(-s['%s_reporter' % state]/rt) * multiplier * reporter
        denom += np.exp(-s[state]/rt) * multiplier + \
            np.exp(-s['%s_reporter' % state]/rt) * multiplier * reporter
    return num/denom


def plot(s, outfile, nsteps=25, reporter=1e-9, title=None):
    """
    make plot with 2D slices in each dimension and one for num vs denom

    Args:
        s (pd.Series): contains fields with free energies of each state
        outfile (str): filename to output to
        nsteps (int): number of increments along each axis
        reporter (float): concentration of reporter
    """
    plt.rc('text')

    m = get_conc_space(s, nsteps, reporter)
    concrange = np.logspace(-15, -3, nsteps)
    ticks = np.arange(0, nsteps, nsteps//4)

    # find slice closest to reporter concentration
    sliceidx = np.argmin(np.abs(np.log10(concrange) - np.log10(reporter)))

    # colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list('earth',
        ['#d7b18e', '#A46C39', '#A48535', '#C2B02E', '#A7A036', '#869B40',
         '#548540', '#3F7F60', '#28726E', '#008299', '#004999', '#002899'])

    # A vs B
    plt.subplot(221)
    plt.imshow(m[:,:,sliceidx].T, cmap=cmap, interpolation='none',
               origin='lower')
    plt.colorbar()
    plt.xlabel('[A]')
    plt.ylabel('[B]')
    plt.xticks(ticks, ['10$^{%d}$' % np.log10(concrange[i]) for i in ticks])
    plt.yticks(ticks, ['10$^{%d}$' % np.log10(concrange[i]) for i in ticks])

    # A vs C
    plt.subplot(222)
    plt.imshow(m[:,sliceidx,:].T, cmap=cmap, interpolation='none',
               origin='lower')
    plt.colorbar()
    plt.xlabel('[A]')
    plt.ylabel('[C]')
    plt.xticks(ticks, ['10$^{%d}$' % np.log10(concrange[i]) for i in ticks])
    plt.yticks(ticks, ['10$^{%d}$' % np.log10(concrange[i]) for i in ticks])

    # B vs C
    plt.subplot(223)
    plt.imshow(m[sliceidx,:,:].T, cmap=cmap, interpolation='none',
               origin='lower')
    plt.colorbar()
    plt.xlabel('[B]')
    plt.ylabel('[C]')
    plt.xticks(ticks, ['10$^{%d}$' % np.log10(concrange[i]) for i in ticks])
    plt.yticks(ticks, ['10$^{%d}$' % np.log10(concrange[i]) for i in ticks])

    m = get_2d_space(s, nsteps, reporter)
    concrange = np.logspace(-30, -6, nsteps)

    # AB vs CC
    plt.subplot(224)
    plt.imshow(m.T, cmap=cmap, interpolation='none', origin='lower')
    plt.colorbar()
    plt.xlabel('[A][B]')
    plt.ylabel('[C]$^2$')
    plt.xticks(ticks, ['10$^{%d}$' % np.log10(concrange[i]) for i in ticks])
    plt.yticks(ticks, ['10$^{%d}$' % np.log10(concrange[i]) for i in ticks])

    if title is not None:
        plt.suptitle(title)
    plt.tight_layout(rect=[0,0,1,.95])
    plt.savefig(outfile, pil_kwargs={"height": 800,"width": 600 })
    plt.close()

def main():
    p = argparse.ArgumentParser()
    p.add_argument('sequence', help='sequence to be simulated')
    p.add_argument('inputfile', help='name of file containing conditions')
    p.add_argument('-o', '--outfile', help='name of output file',
                   default='concentration_heatmap.png')
    p.add_argument('-t', '--title', help='title at top of plot')
    p.add_argument('-n', '--nsteps', default=25, type=int,
                   help='number of increments along each axis')
    args = p.parse_args()

    inputs, reporter, complexes = parse_states(args.inputfile)
    s = get_pfs(args.sequence, inputs, reporter, complexes)
    plot(s, args.outfile, nsteps=args.nsteps, title=args.title)


if __name__ == '__main__':
    main()
