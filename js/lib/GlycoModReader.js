/**
 * @fileOverview    Classes for reading GlycoMod data
 */

if ( typeof MASCP == 'undefined' || typeof MASCP.Service == 'undefined' ) {
    throw "MASCP.Service is not defined, required class";
}

/** Default class constructor
 *  @class      Service class that will retrieve data from GlycoMod for a given AGI.
 *  @param      {String} agi            Agi to look up
 *  @param      {String} endpointURL    Endpoint URL for this service
 *  @extends    MASCP.Service
 */
MASCP.GlycoModReader = MASCP.buildService(function(data) {
                        this._raw_data = data;
                        return this;
                    });

MASCP.GlycoModReader.SERVICE_URL = '/?';

MASCP.GlycoModReader.prototype.requestData = function()
{
    var agi = this.agi;

    return {
        type: "GET",
        dataType: "json",
        data: { 'agi'       : agi,
                'service'   : 'glycomod'
        }
    };
};

/**
 *  @class   Container class for results from the GlycoMod service
 *  @extends MASCP.Service.Result
 */
// We need this line for the JsDoc to pick up this class
MASCP.GlycoModReader.Result = MASCP.GlycoModReader.Result;

/** Retrieve the peptides for this particular entry from the GlycoMod service
 *  @returns Array of peptide strings
 *  @type [String]
 */
MASCP.GlycoModReader.Result.prototype.getPeptides = function()
{
    var content = null;

    if (this._peptides) {
        return this._peptides;
    }

    if (! this._raw_data || ! this._raw_data.peptides) {
        return [];
    }


    this._peptides= this._raw_data.peptides;

    return this._peptides;
};

MASCP.GlycoModReader.Result.prototype._cleanSequence = function(sequence)
{
    return sequence.replace(/[^A-Z]/g,'');
};

MASCP.GlycoModReader.prototype.setupSequenceRenderer = function(sequenceRenderer)
{
    var reader = this;

    var css_block = '.active .overlay { background: #ff00ff; } .active a { color: #000000; text-decoration: none !important; }  :indeterminate { background: #ff0000; } .tracks .active { background: #0000ff; } .inactive a { text-decoration: none; } .inactive { display: none; }';

    this.bind('resultReceived', function() {
        var peps = this.result.getPeptides();

        var overlay_name = 'glycomod_experimental';
        var icons = [];

        if (peps.length > 0) {
            MASCP.registerLayer(overlay_name,{ 'fullname' : 'Glycosylation (mod)', 'color' : '#ff00ff', 'css' : css_block, 'hover_peptides' : true });

            MASCP.registerGroup('glycomod_peptides', {'fullname' : 'Glycosilation Publication ', 'hide_member_controllers' : true, 'hide_group_controller' : true, 'color' : '#ff00ff' });
            if (sequenceRenderer.createGroupController) {
                sequenceRenderer.createGroupController('glycomod_experimental','glycomod_peptides');
            }
        }

        for (var j = 0; j < peps.length; j++ ) {
            var pep = peps[j];

            if (pep.length === 0) {
                continue;
            }
            var layer_name = 'glycomod_spectrum_';
            MASCP.registerLayer(layer_name, { 'fullname': 'Spectrum ', 'group' : 'glycomod_peptides', 'color' : '#ff00ff', 'css' : css_block });
            var peptide = pep.sequence;
            var peptide_bits = sequenceRenderer.getAminoAcidsByPeptide(peptide);
            if (peptide_bits.length === 0){
	        continue;
            }
            peptide_bits.addToLayer(layer_name);
            icons.push(peptide_bits.addToLayer('glycomod_experimental'));

            for (var k = 0; k < pep.de_index.length; k++ ) {
                var bit = peptide_bits[pep.de_index[k] - 1];
    	       	icons = icons.concat(peptide_bits[pep.de_index[k] - 1].addToLayer('glycomod_experimental', { 'content' : 'GLY', 'height' : 20, 'offset' : -2.5, 'popup': { 'Glycosylation': bit.amino_acid+' '+bit._index } }));
        		peptide_bits[pep.de_index[k] - 1].addToLayer(layer_name, {'content' : 'GLY', 'height' : 20, 'offset' : -2.5 });
            }
        }
        jQuery(sequenceRenderer).trigger('resultsRendered',[reader]);
    });
    return this;
};

/** Retrieve an array of positions that phosphorylation has been experimentally verified to occur upon
 *  @returns {Array}    Phosphorylation positions upon the full protein
 */
MASCP.GlycoModReader.Result.prototype.getAllExperimentalPositions = function()
{
    var specs = this.getPeptides();
    var results = [];
    var seen = {};
    specs.forEach(function(spec) {
        var peps = spec.peptides;
        peps.forEach(function(pep) {
            pep.de_index.forEach(function(pos) {
                if ( ! seen[pos] ) {
                    results.push(pos);
                    seen[pos] = true;
                }
            });
        });
    });
    return results;
}

MASCP.GlycoModReader.Result.prototype.render = function()
{
};

