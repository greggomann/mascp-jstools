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
MASCP.Modhunter.prototype.loadSequence = function(modObject, inputSequence) {
    
    // Initialize the sequence for the Modhunter
    // Whole protein sequence is needed later to find peptide indeces
    modObject['whole_sequence'] = inputSequence;
    modObject['predicted_total'] = 0;
    modObject['peptide_total'] = 0;

    for (var i = 0; i < inputSequence.length; i++) {

        // Residue index in the protein sequence is also the residue sequence here
        modObject.sequence[i] = {};

        // Coverage properties are used by modhunter to compute scores for each residue
        modObject.sequence[i]['gator_coverage'] = 0;
        modObject.sequence[i]['predicted_coverage'] = 0;
        modObject.sequence[i]['reader_coverage'] = 0;
        modObject.sequence[i]['score'] = 0;
    }
};


// countCoverage method counts the number of peptides covering each residue
// This method is configured to count peptides from the readers contained in readerList object
MASCP.Modhunter.prototype.countCoverage = function(modObject, reader) {

    // readerList contains the names of readers to be included in the count
    //      along with an array containing four parameters, for example:
    //          { 'MASCP.ReaderName' : [count_fieldname, count_type_flag, norm_factor] }
    //          count_fieldname == name of field in peptide object containing # of experiments/spectral count
    //          count_type_flag == 'length' if peptide count is given by len(count_fieldname)
    //          count_type_flag == 'value' if peptide count is given by the value of count_fieldname
    //          norm_factor == normalization factor for peptide count in each database
    var readerList = { 'MASCP.PpdbReader': ['experiments', 'length', 0.00732903204172],
                    'MASCP.AtPeptideReader': ['tissues', 'length', 0.00277509729109],
                    'MASCP.Pep2ProReader': ['qty_spectra', 'value', 0.00229551940475],
                    'MASCP.AtChloroReader': ['', '', 0.0570994709686],
                    'MASCP.GelMapReader': ['', '', 1.0],
                    'MASCP.ProteotypicReader': ['pvalue', 'value', 1.0],
                    'MASCP.PubmedReader': ['', '', 0.321208781115]  };

    if (typeof readerList[reader.toString()] === 'undefined') {
        return;
    }
    
    // Function to convert a list of peptides into a list of sequence indeces
    var convertToIndeces = function(pepSeq) {
        var results = [];
        var hitIndex = modObject.whole_sequence.indexOf(pepSeq);
        if (hitIndex >= 0) {
            results = [hitIndex, (hitIndex+pepSeq.length)];
        } else {
            results = [0,0];
        }

        return results;
    };
    
    // Populate list of peptides for which to count coverage
    var getPeps = [];
    if (reader.result) {
        getPeps = reader.result.getPeptides();
    }
    
    // rdrCoverageList keeps track of indeces whose reader_coverage property have
    //      already been incremented for the current reader
    var rdrCoverageList = [];
    
    var incCoverage = function(propName, value, leftIndex, rightIndex) {
        for (var i = leftIndex; i < rightIndex; i++) {
            modObject.sequence[p][propName] += value;
        }
    };

    if (getPeps.length > 0) {
        // For each peptide, increment coverage properties for the residues within the peptide
        for (var pep in getPeps) {
            var thisIdx = convertToIndeces(getPeps[pep].sequence);
            for (var p = thisIdx[0]; p < thisIdx[1]; p++) {
                // Check for predicted peptide vs. experimental peptide
                if (reader.toString() == 'MASCP.ProteotypicReader') {
                    modObject.sequence[p].predicted_coverage += 1;
                    modObject.predicted_total += 1;
                } else {
                    // Retrieve spectral count for experimental peptides, if available
                    var pepCount = 1;
                    if (readerList[reader.toString()][0] == 'length') {
                        pepCount = getPeps[pep][readerList[reader.toString()][0]].length;
                    }
                    else if (readerList[reader.toString()][0] == 'value') {
                        pepCount = parseInt(getPeps[pep][readerList[reader.toString()][0]]);
                    }
                    modObject.peptide_total += pepCount;
                    modObject.sequence[p].gator_coverage += pepCount;
                    // If modObject residue's reader_coverage hasn't yet been incremented, do so now
                    if (rdrCoverageList.indexOf(p) < 0) {
                        modObject.sequence[p].reader_coverage += 1;
                        rdrCoverageList.push(p);
                    }
                }
            }
        }
    }
};


// calcScores method calculates Modhunter scores for each residue
MASCP.Modhunter.prototype.calcScores = function() {
    
    var seqLength = this.whole_sequence.length;

    // Function to convert a sequence index to a percentage
    // Returned value is a floating point between 0 and 1
    var toPercent = function(index) {
        return Math.min((index / seqLength), 1);
    };

    // Set maximum confidence of modification suggestions based on total # of peptides in GATOR
    // var maxRating = Math.min((this.peptide_total / 15), 1);
    
    // Iterate through amino acids and compute modhunter rating from 0-100 for each
    for (var q = 0; q < seqLength; q++) {
        var thisRating = 100;
        var predictedBool = 0;
        var readerDec = (this.sequence[q].reader_coverage < 2) ? 0.5 : 1;
        // var gatorDec = Math.max(1 - (this.sequence[q].gator_coverage / 4), 0);
        if (this.sequence[q].predicted_coverage > 0) { predictedBool = 1; }
        thisRating = thisRating * predictedBool;
        // thisRating = thisRating * maxRating;
        this.sequence[q].score = thisRating;
    }

};
