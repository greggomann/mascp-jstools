/*
Functions and data structures to evaluate the Gator coverage of residues in a protein
and compute Modhunter scores for each residue

Author: Greg Mann
JBEI, 12/8/24
*/

MASCP.Modhunter = function() {
    // sequence property holds residue indeces and associated variables
    this.sequence = {};
    // confirmed_mods is an object literal containing experimentally-confirmed modifications from the Gator
    this.confirmed_mods = [];
};


/**
 *  Binds a handler to one or more events. Returns a reference to self, so this method
 *  can be chained.
 *
 *  @param  {String}    type        Event type to bind
 *  @param  {Function}  function    Handler to execute on event
 */

MASCP.Modhunter.prototype.bind = function(type,func)
{
    bean.add(this,type,func);
    return this;
};


// loadSequence method initializes Modhunter objects for a given protein sequence
MASCP.Modhunter.prototype.loadSequence = function(modObject, inputSequence) {

    // Initialize the sequence for the Modhunter
    // Whole protein sequence is needed later to find peptide indeces
    modObject['whole_sequence'] = (inputSequence) ? inputSequence : '';
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

    bean.fire(this,'sequenceLoaded');
    return this;
};


// countCoverage method counts the number of peptides covering each residue
// This method is configured to count peptides from the readers contained in readerList object
MASCP.Modhunter.prototype.countCoverage = function(modObject, reader) {
    var self = this;
    // readerList contains the names of readers to be included in the count
    //      along with an array containing four parameters, for example:
    //          { 'MASCP.ReaderName' : [count_fieldname, count_type_flag, norm_factor] }
    //          count_fieldname == name of field in peptide object containing # of experiments/spectral count
    //          count_type_flag == 'length' if peptide count is given by len(count_fieldname)
    //          count_type_flag == 'value' if peptide count is given by the value of count_fieldname
    var readerList = { 'MASCP.PpdbReader': ['experiments', 'length'],
                    'MASCP.AtPeptideReader': ['tissues', 'length'],
                    'MASCP.Pep2ProReader': ['qty_spectra', 'value'],
                    'MASCP.AtChloroReader': ['', ''],
                    'MASCP.GelMapReader': ['', ''],
                    'MASCP.ProteotypicReader': ['pvalue', 'value'],
                    'MASCP.PubmedReader': ['', ''],
                    'MASCP.SnpReader': ['', ''],
                    'MASCP.OrthologyReader': ['', ''] };

    var getIndex = function(pepSeq) {
        return self.whole_sequence.indexOf(pepSeq);
    }

    // modRdrList contains the names of readers that provide experimentally
    //      confirmed modification sites
    //          { 'MASCP.ReaderName': 
    var modRdrList = { 'MASCP.RnaEditReader': function() {
            var mods = [];
            var accessions = this.getAccessions();
            while (accessions.length > 0) {
                var acc = accessions.shift();
                var edits = this.getSnp(acc);
                for (var i = 0; i < edits.length; i++) {
                    mods.push([edits[i][0], 'RNA Edit ('+edits[i][1]+' to '+edits[i][2]+')']);
                }
            }
            return mods;
        }, 'MASCP.GlycoModReader': function() {
            var mods = [];
            var peps = this.getPeptides();
            for (var i = 0; i < peps.length; i++) {
                var pepIdx = getIndex(peps[i].sequence);
                for (var j = 0; j < peps[i].de_index.length; j++) {
                    mods.push([pepIdx+peps[i].de_index[j]-1, 'Glycosylation']);
                }
            }
            return mods;
        }, 'MASCP.PhosphatReader': function() {
            var mods = [];
            var sites = this.getAllExperimentalPositions();
            for (var i = 0; i < sites.length; i++) {
                mods.push([sites[i]-1, 'Phosphorylation']);
            }
            return mods;
        }, 'MASCP.UbiquitinReader': function() {
            var mods = [];
            var peps = this.getPeptides();
            for (var i = 0; i < peps.length; i++) {
                var pepIdx = self.whole_sequence.indexOf(peps[i].sequence);
                for (var j = 0; j < peps[i].positions.length; j++) {
                    mods.push([pepIdx+peps[i].positions[j]-1, 'Ubiquitination']);
                }
            }
            return mods;
        }, 'MASCP.RippdbReader': function() {
            var mods = [];
            var spectra = this.getSpectra();
            for (var i = 0; i < spectra.length; i++) {
                for (var j = 0; j < spectra[i].peptides.length; j++) {
                    var pepIdx = self.whole_sequence.indexOf(spectra[i].peptides[j].sequence);
                    for (var k = 0; k < spectra[i].peptides[j].positions.length; k++) {
                        mods.push([pepIdx+spectra[i].peptides[j].positions[k]-1, 'Phosphorylation']);
                    }
                }
            }
            return mods;
        }
    };

    var readerName = reader.toString();
    var isRdr = !(typeof readerList[readerName] === 'undefined');
    var isModRdr = !(typeof modRdrList[readerName] === 'undefined');
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
    
    // Populate list of peptides, modification sites, or nsSNPs for which to count coverage
    var getPeps = [];
    if (reader.result) {
        // If this is a normal Gator reader with standard shotgun peptides, do the following
        if (isRdr) {
            // If this is for the SnpReader, use getSnp() method instead of getPeptides()
            if (readerName == 'MASCP.SnpReader') {
                var accList = reader.result.getAccessions();
                for (var accIdx in accList) {
                    getPeps.push.apply(getPeps, reader.result.getSnp(accList[accIdx]));
                }
            } else if (readerName == 'MASCP.OrthologyReader') {
                getPeps = reader.result.getConservation();
            } else {
                getPeps = reader.result.getPeptides();
            }
        // If this is a reader containing confirmed modification data, do the following
        } else if (isModRdr) {
            var modResults = modRdrList[readerName].call(reader.result);
            this.confirmed_mods = this.confirmed_mods.concat(modRdrList[readerName].call(reader.result));
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
            switch (readerName) {
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
                    if (modObject.sequence[getPeps[pep][0]-1]) {
                        modObject.sequence[getPeps[pep][0]-1].snp_coverage += 1;
                    }
                    break;
                case 'MASCP.OrthologyReader':
                    if (modObject.sequence[pep]) {
                        modObject.sequence[pep].conservation = getPeps[pep];
                    }
                default:
                    var thisIdx = convertToIndeces(getPeps[pep].sequence);
                    var firstLoop = true;
                    for (var p = thisIdx[0]; p < thisIdx[1]; p++) {
                        // Retrieve spectral count for experimental peptides, if available
                        var pepCount = 1;
                        if (readerList[readerName][1] == 'length') {
                            pepCount = getPeps[pep][readerList[readerName][0]].length;
                        }
                        else if (readerList[readerName][1] == 'value') {
                            pepCount = parseInt(getPeps[pep][readerList[readerName][0]]);
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

    bean.fire(this,'coverageCounted');
    return this;
};


// calcScores method calculates two different scores:
// Modhunter scores for each residue, and the abundance score for the protein
MASCP.Modhunter.prototype.calcScores = function() {
    // binMap contains log-scale bin values for the protein abundance score
    var binMap = [1.0, 1.11073537957, 1.23373308343, 1.37035098472, 1.52209732116, 1.69064734577, 1.87786182132, 2.08580756289, 2.31678025509, 2.57332979602, 2.85828844775, 3.17480210394, 3.52636501998, 3.91685838898, 4.35059318942, 4.83235777762, 5.36747075035, 5.96183966124, 6.62202623908, 7.3553188282, 8.16981285052, 9.07450017756, 10.0793683992, 11.1955110847, 12.4352502542, 13.8122724111, 15.3417796394, 17.040657431, 18.9276610998, 21.0236228361, 23.3516816909, 25.9375390266, 28.8097422559, 32.0, 35.5435321463, 39.4794586699, 43.851231511, 48.7071142772, 54.1007150645, 60.0915782823, 66.7458420126, 74.1369681627, 82.3465534726, 91.4652303279, 101.593667326, 112.843680639, 125.339468447, 139.218982061, 154.635448884, 171.759064011, 190.77886916, 211.904839651, 235.370202503, 261.434011217, 290.384005682, 322.539788773, 358.25635471, 397.928008133, 441.992717157, 490.936948459, 545.301037793, 605.685155195, 672.755930757, 747.253814109, 830.001248851, 921.911752189, 1024.0, 1137.39302868, 1263.34267744, 1403.23940835, 1558.62765687, 1731.22288206, 1922.93050504, 2135.8669444, 2372.38298121, 2635.08971112, 2926.88737049, 3250.99735443, 3610.99778046, 4010.86299032, 4455.00742597, 4948.33436428, 5496.29004836, 6104.92381311, 6780.95486882, 7531.84648008, 8365.88835894, 9292.28818183, 10321.2732407, 11464.2033507, 12733.6962603, 14143.766949, 15709.9823507, 17449.6332094, 19381.9249662, 21528.1897842, 23912.1220515, 26560.0399632, 29501.17607, 32768.0];

    var seqLength = this.whole_sequence.length;

    // totalCoverage contains the sum of gator_coverage over all residues
    var totalCoverage = 0;
    for (var j = 0; j < seqLength; j++) {
        totalCoverage += this.sequence[j].gator_coverage;
    }

    // Set protein abundance score
    var abScore = (this.peptide_total / this.tryptic_total) * 100;
    
    // Find appropriate bin for this abundance score, giving us the final abundance score
    var i = 0;
    while (binMap[i] < abScore && i < 99) {
        i++;
    }
    this.abundance_score = i;
    jQuery('#abundance').text(i);

    // Look for signs of C or N-terminal processing
    if (this.abundance_score >= 50) {
        var covCount = 0;
        var termLength = (seqLength > 200) ? seqLength * 0.2 : 40;
        for (var k = 0; covCount < (totalCoverage*0.05); k++) {
            covCount += this.sequence[k].gator_coverage;
        }
        if (k > parseInt(seqLength*0.08) && k < termLength) {
            // n_terminal property contains the residue index at end of processing region
            this.n_terminal = k-1;
        }
        covCount = 0;
        for (var k = seqLength-1; covCount < (totalCoverage*0.05); k--) {
            covCount += this.sequence[k].gator_coverage;
        }
        if (seqLength-k > parseInt(seqLength*0.08) && k < termLength) {
            // c_terminal property contains residue index at beginning of processing region
            this.c_terminal = k+1;
        }
    }

    // Iterate through amino acids and compute modhunter rating from 0-100 for each
    for (var q = 0; q < seqLength; q++) {
        // gatScore is based on # of peptides in the Gator that cover this residue
        var gatScore = Math.min(this.sequence[q].gator_coverage / 10, 1);
        // predScore is based on # of predicted peptides that cover this residue
        var predScore = (this.sequence[q].predicted_coverage > 0) ? 1 : 0;
        // abScale scales the ModHunter score based on protein abundance score
        var abScale = Math.min(Math.max(this.abundance_score - 20, 0) / 50, 1);
        // gapScale adds to the mod score in gaps when protein abundance is high
        var gapScale = (Math.max(this.abundance_score - 70, 0) / 30) * Math.max((4 - this.sequence[q].gator_coverage) / 4, 0) * 0.7;
        // conScale is based on the orthology track's conservation value
        // var conScale = (this.sequence[q].conservation) ? this.sequence[q].conservation : 1;
        // modScore is the ModHunter score
        var modScore = Math.round(Math.min(Math.max(predScore - gatScore + gapScale, 0), 1) * abScale * 100);
        if (this.n_terminal && q <= this.n_terminal) {
            modScore = Math.min(modScore+60, 100);
        } else if (this.c_terminal && q >= this.c_terminal) {
            modScore = Math.min(modScore+60, 100);
        }
        this.sequence[q].score = modScore;
    }

    bean.fire(this,'scoresCalculated');
    return this;
};
