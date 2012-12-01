'''
Created on 30 Oct 2012

@author: musselle
'''
 while not input_checked_flag:
                try:
                    ans = raw_input('Proceed? [y/n]: ')
                    if ans == 'n':
                        sys.exit('Quit by user.')
                    elif ans == 'y':
                        idxFileName = filename
                        input_checked_flag = 0
                        break
                    else:
                        raise ValueError('ans')
                except ValueError:
                    print 'Invalid input. Try again...'
            
            try:
                yield SeqIO.index_db(idxFileName, format='fastq')
            except IOError as e:
                print e
                print 'Invalid file name, or {0} may not exist'.format(filename)
        