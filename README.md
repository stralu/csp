# CSP 1.1
Here you can find a fork of the original CSP v1.1 plugin for EEGLAB @P. Sajda that is available at: [EEGLAB Extensions](https://sccn.ucsd.edu/wiki/EEGLAB_Extensions).

In the original version it can happen that too few trials are used for CSP calculation or that the trial-index exceeds the size of the data matrix. This occurs, if ALLEEG consists of more than two EEG structures and CSP is not calculated between the first two.

This is fixed here.
