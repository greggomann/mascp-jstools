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
    modObject['abundance_score'] = 0;

    var trypticCount = 1;

    for (var i = 0; i < inputSequence.length; i++) {
        // Compute in silico tryptic digest
        if (i != inputSequence.length-1 && inputSequence[i] in {'K':'','R':''} && inputSequence[i+1] != 'P') {
            trypticCount++;
        }

        // Residue index in the protein sequence is also the residue sequence here
        modObject.sequence[i] = {};

        // Coverage properties are used by modhunter to compute scores for each residue
        modObject.sequence[i]['gator_coverage'] = 0;
        modObject.sequence[i]['predicted_coverage'] = 0;
        modObject.sequence[i]['reader_coverage'] = 0;
        modObject.sequence[i]['snp_coverage'] = 0;
        modObject.sequence[i]['score'] = 0;
    }

    modObject['tryptic_total'] = trypticCount;
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
                    'MASCP.PubmedReader': ['', '', 0.321208781115],
                    'MASCP.SnpReader': ['', '', 1.0] };

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
        // If this is for the SnpReader, use getSnp() method instead of getPeptides()
        if (reader.toString() == 'MASCP.SnpReader') {
            var accList = reader.result.getAccessions();
            for (var accIdx in accList) {
                getPeps.push.apply(getPeps, reader.result.getSnp(accList[accIdx]));
            }
        } else {
            getPeps = reader.result.getPeptides();
        }
    }
    
    // rdrCoverageList keeps track of indeces whose reader_coverage property have
    //      already been incremented for the current reader
    var rdrCoverageList = [];

    if (getPeps.length > 0) {
        // For each peptide, increment coverage properties for the residues covered by the peptide
        for (var pep in getPeps) {
            // Special cases for ProteotypicReader or SnpReader
            // Otherwise, run default routine for generic reader
            switch (reader.toString()) {
                case 'MASCP.ProteotypicReader':
                    var thisIdx = convertToIndeces(getPeps[pep].sequence);
                    var firstLoop = true;
                    for (var p = thisIdx[0]; p < thisIdx[1]; p++) {
                        modObject.sequence[p].predicted_coverage += 1;
                        if (firstLoop == true) {
                            modObject.predicted_total += 1;
                            firstLoop = false;
                        }
                    }
                    break;
                case 'MASCP.SnpReader':
                    modObject.sequence[getPeps[pep][0]-1].snp_coverage += 1;
                    break;
                default:
                    var thisIdx = convertToIndeces(getPeps[pep].sequence);
                    var firstLoop = true;
                    for (var p = thisIdx[0]; p < thisIdx[1]; p++) {
                        // Retrieve spectral count for experimental peptides, if available
                        var pepCount = 1;
                        if (readerList[reader.toString()][1] == 'length') {
                            pepCount = getPeps[pep][readerList[reader.toString()][0]].length;
                        }
                        else if (readerList[reader.toString()][1] == 'value') {
                            pepCount = parseInt(getPeps[pep][readerList[reader.toString()][0]]);
                        }
                        if (firstLoop == true) {
                            modObject.peptide_total += pepCount;
                            firstLoop = false;
                        }
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


// calcScores method calculates Modhunter scores for each residue and the abundance score for the protein
MASCP.Modhunter.prototype.calcScores = function() {
    // binMap contains log-scale bin values for the protein abundance score
    var binMap = [1.0, 1.11073537957, 1.23373308343, 1.37035098472, 1.52209732116, 1.69064734577, 1.87786182132, 2.08580756289, 2.31678025509, 2.57332979602, 2.85828844775, 3.17480210394, 3.52636501998, 3.91685838898, 4.35059318942, 4.83235777762, 5.36747075035, 5.96183966124, 6.62202623908, 7.3553188282, 8.16981285052, 9.07450017756, 10.0793683992, 11.1955110847, 12.4352502542, 13.8122724111, 15.3417796394, 17.040657431, 18.9276610998, 21.0236228361, 23.3516816909, 25.9375390266, 28.8097422559, 32.0, 35.5435321463, 39.4794586699, 43.851231511, 48.7071142772, 54.1007150645, 60.0915782823, 66.7458420126, 74.1369681627, 82.3465534726, 91.4652303279, 101.593667326, 112.843680639, 125.339468447, 139.218982061, 154.635448884, 171.759064011, 190.77886916, 211.904839651, 235.370202503, 261.434011217, 290.384005682, 322.539788773, 358.25635471, 397.928008133, 441.992717157, 490.936948459, 545.301037793, 605.685155195, 672.755930757, 747.253814109, 830.001248851, 921.911752189, 1024.0, 1137.39302868, 1263.34267744, 1403.23940835, 1558.62765687, 1731.22288206, 1922.93050504, 2135.8669444, 2372.38298121, 2635.08971112, 2926.88737049, 3250.99735443, 3610.99778046, 4010.86299032, 4455.00742597, 4948.33436428, 5496.29004836, 6104.92381311, 6780.95486882, 7531.84648008, 8365.88835894, 9292.28818183, 10321.2732407, 11464.2033507, 12733.6962603, 14143.766949, 15709.9823507, 17449.6332094, 19381.9249662, 21528.1897842, 23912.1220515, 26560.0399632, 29501.17607, 32768.0];

    var seqLength = this.whole_sequence.length;

    // Set protein abundance score
    var abScore = (this.peptide_total / this.tryptic_total) * 100;
    
    // Find appropriate bin for this abundance score, giving us the final abundance score
    var i = 0;
    while (binMap[i] < abScore && i < 99) {
        i++;
    }
    this.abundance_score = i;
    jQuery('#abundance').text(i);

    // Iterate through amino acids and compute modhunter rating from 0-100 for each
    for (var q = 0; q < seqLength; q++) {
        // gatScore is based on # of peptides in the Gator that cover this residue
        var gatScore = Math.max(6 - this.sequence[q].gator_coverage, 0) / 6;
        // predScore is based on # of predicted peptides that cover this residue
        var predScore = (this.sequence[q].predicted_coverage > 0) ? 1 : 0;
        // snpScore is based on # of nsSNPs
        var snpScore = Math.max(1-(this.sequence[q].snp_coverage / 50), 0);
        // abScale scales the ModHunter score based on protein abundance score
        var abScale = Math.min(this.abundance_score / 60, 1);
        // modScore is the ModHunter score
        var modScore = Math.round(Math.min(Math.round(((gatScore * 60) + (predScore * 40)) * snpScore), 100) * abScale);
        this.sequence[q].score = modScore;
    }
};
