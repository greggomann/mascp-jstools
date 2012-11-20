/*
Functions and data structures to evaluate the Gator coverage of residues in a protein
and compute Modhunter scores for each residue

Author: Greg Mann
JBEI, 12/8/24
*/

MASCP.Modhunter = function() {
    // sequence property holds residue indeces and associated variables
    this.sequence = {};
};


// loadSequence method initializes Modhunter objects for a given protein sequence
MASCP.Modhunter.prototype.loadSequence = function(inputSequence) {
    
    // Initialize the sequence for the Modhunter
    for (var i = 0; i < inputSequence.length; i++) {

        // Whole protein sequence is needed later to find peptide indeces
        this.sequence['whole_sequence'] = inputSequence;

        // Residue index in the protein sequence is also the residue sequence here
        this.sequence[i] = {};

        // Coverage properties are used by modhunter to compute scores for each residue
        this.sequence[i]['gator_coverage'] = 0;
        this.sequence[i]['predicted_coverage'] = 0;
        this.sequence[i]['reader_coverage'] = 0;
        this.sequence[i]['modhunter_score'] = 0;
    }
};


// countCoverage method counts the number of peptides covering each residue
// This method is configured to count peptides from the readers contained in readerList object
MASCP.Modhunter.prototype.countCoverage = function(reader) {
    var readerList = { 'MASCP.PpdbReader': 'sequence',
                    'MASCP.AtPeptideReader': 'sequence',
                    'MASCP.Pep2ProReader': '',
                    'MASCP.AtChloroReader': 'sequence',
                    'MASCP.GelMapReader': '',
                    'MASCP.ProteotypicReader': 'sequence' };

    // Function to convert a list of peptides into a list of sequence indeces
    var convertToIndeces = function(peptideList) {
        var results = [];

        for (var i = 0; i < peptideList.length; i++) {
            var hitIndex = this.sequence['whole_sequence'].indexOf(peptideList[i]);
            if (hitIndex >= 0) {
                results.push([hitIndex, (hitIndex+peptideList[i].length)]);
            } else {
                results.push([0,0]);
            }
        }

        return results;
    };

    // Populate list of peptides to count coverage for
    pepList = [];
    if (reader.result && readerList[reader.toString()]) {
        getPeps = reader.result.getPeptides();
        if (getPeps) {
            // Refer to readerList to see if the list of peptides is given directly or
            //     if a subfield must be used
            if (readerList[reader.toString()].pepField == '') {
                pepList = pepList.concat(getPeps);
            } else {
                for (j = 0; j < getPeps.length; j++) {
                    pepList.push(getPeps[j][readerList[reader.toString()].pepField]);
                }
            }
        }
    }

    // rdrCoverageList keeps track of indeces whose reader_coverage property have
    //      already been incremented for this reader
    var rdrCoverageList = [];
    if (pepList.length > 0) {
        pepList = convertToIndeces(pepList);
        // For each peptide, increment coverage properties for the residues within the peptide
        for (var o = 0; o < pepList.length; o++) {
            for (var p = pepList[o][0]; p < pepList[o][1]; p++) {
                if (reader.toString() == 'MASCP.ProteotypicReader') {
                    this.sequence[p].predicted_coverage += 1;
                } else {
                    this.sequence[p].gator_coverage += 1;
                    // If this residue's reader_coverage hasn't yet been incremented, do so now
                    if (rdrCoverageList.indexOf(p) < 0) {
                        this.sequence[p].reader_coverage += 1;
                        rdrCoverageList.push(p);
                    }
                }
            }
        }
    }

};


// calcScores method calculates Modhunter scores for each residue
MASCP.Modhunter.prototype.calcScores = function() {

    var sequence = this.sequence;
    var seqLength = sequence.length;

    // Function to convert a sequence index to a percentage
    // Returned value is a floating point between 0 and 1
    var toPercent = function(index) {
        return Math.min((index / seqLength), 1);
    };

    // Set maximum confidence of modification suggestions based on total # of peptides in GATOR
    var maxRating = numPeps / 40;
    maxRating = (maxRating > 1) ? 1 : maxRating;

    // Iterate through amino acids and compute modhunter rating from 0-100 for each
    for (var q = 0; q < this._sequence_els.length; q++) {
        stopObject[q] = 0;
        var thisRating = 0;
        thisRating += (4 - this._sequence_els[q]._gator_coverage) * 9;
        thisRating += this._sequence_els[q]._predicted_coverage * 22;
        thisRating -= this._sequence_els[q]._reader_coverage * 20;
        thisRating = Math.max(Math.min(thisRating, 100), 0);
        thisRating = thisRating * maxRating;
        stopObject[q] = thisRating;
    }

};
