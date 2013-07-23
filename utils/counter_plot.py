'''
Created on Jul 23, 2013

@author: Chris Musselle
'''

from collections import Counter

def plot_counters(counters, labels=None, log='xy', xlab="", ylab="", title="", **kwargs):
    ''' Construct a series of scatter plots from a list of Counter Dictionaries '''
    
    import matplotlib.pyplot as plt
    
    # Default plot options
    if 'ms' not in kwargs:   # set marker size
        kwargs['ms'] = 4.0
    if 'marker' not in kwargs: # set marker type
        kwargs['marker'] = '.'
    if 'mew' not in kwargs: # set marker edge weight 
        kwargs['mew'] = 0.0 
    
    if type(counters) is not list and type(counters) is not tuple:
        counters = [counters]
    
    if labels is not None:
        assert len(labels) == len(counters), "Number of labels must match number of counters."
    
    for i in range(len(counters)):
        
        data_xs = [int(k) for k in counters[i].keys()]
        data_ys = counters[i].values()

        if labels:
            plt.plot(data_xs, data_ys, label=labels[i], ls='', **kwargs)
        else:
            plt.plot(data_xs, data_ys, label='Counter-'+str(i), ls='', **kwargs)
        
        ax = plt.gca()
        if 'x' in log:
            ax.set_xscale('log')
        if 'y' in log:
            ax.set_yscale('log')
        
        plt.title(title)
        plt.xlabel(xlab)
        plt.ylabel(ylab)

    plt.legend(numpoints=1, markerscale=8)
    plt.show()


if __name__ == '__main__' : 
    
    filename = ''
    
    f = open(filename, 'rb')
    counter = Counter() 

    for line in f:
        #
        # insert code to get correct bin for line 
        #
        
        counter[bin] +=1 
        
    
    plot_counters(counter, log='', xlab='Bins', ylab='Frequency')


    
     



