/** @fileOverview   Classes for reading data from the NTerm data
 */
if ( typeof MASCP == 'undefined' || typeof MASCP.Service == 'undefined' ) {
    throw "MASCP.Service is not defined, required class";
}

/** Default class constructor
 *  @class      Service class that will retrieve data from the NTerm data for a given AGI.
 *  @param      {String} agi            Agi to look up
 *  @param      {String} endpointURL    Endpoint URL for this service
 *  @extends    MASCP.Service
 */
MASCP.NTermReader = MASCP.buildService(function(data) {
                        this._raw_data = data;
                        return this;
                    });

MASCP.NTermReader.SERVICE_URL = '?';

MASCP.NTermReader.prototype.requestData = function()
{
    var agi = this.agi;

    return {
        type: "GET",
        dataType: "json",
        data: { 'agi'       : agi,
                'service'   : 'nterm'
        }
    };
};

/**
 *  @class   Container class for results from the NTerm service
 *  @extends MASCP.Service.Result
 */
// We need this line for the JsDoc to pick up this class
MASCP.NTermReader.Result = MASCP.NTermReader.Result;

/** Retrieve the peptides for this particular entry from the NTerm service
 *  @returns Array of peptide strings
 *  @type [String]
 */
MASCP.NTermReader.Result.prototype.getPeptides = function()
{
    var content = null;
    if (! this._raw_data || ! this._raw_data.data  || ! this._raw_data.data.peptides ) {
        return [];
    }

    return this._raw_data.data.peptides;
};

MASCP.NTermReader.Result.prototype._cleanSequence = function(sequence)
{
    return sequence.replace(/[^A-Z]/g,'');
};

MASCP.NTermReader.prototype.setupSequenceRenderer = function(sequenceRenderer)
{
    var reader = this;

    var css_block = '.active .overlay { background: #666666; } .active a { color: #000000; text-decoration: none !important; }  :indeterminate { background: #ff0000; } .tracks .active { background: #0000ff; } .inactive a { text-decoration: none; } .inactive { display: none; }';

    this.bind('resultReceived', function() {
        var peps = this.result.getPeptides();

        var overlay_name = 'nterm_experimental';
        var group_name = 'nterm_peptides';
        var icons = [];

        if (peps.length > 0) {
            MASCP.registerLayer(overlay_name,{ 'fullname' : 'N-term Process (mod)', 'color' : '#666666', 'css' : css_block });

            MASCP.registerGroup(group_name, {'fullname' : 'N-term Process (mod)', 'hide_member_controllers' : true, 'hide_group_controller' : true, 'color' : '#666666' });
            if (sequenceRenderer.createGroupController) {
                sequenceRenderer.createGroupController(overlay_name,group_name);
            }

            jQuery(MASCP.getGroup(group_name)).bind('visibilityChange',function(e,rend,vis) {
                if (rend != sequenceRenderer) {
                    return;
                }
                icons.forEach(function(el) {
                    el.style.display = vis ? 'none' : 'inline';
                });
            });


        }

        for (var i = 0; i < peps.length; i++) {
            var layer_name = 'nterm_peptide_'+i;
            MASCP.registerLayer(layer_name, { 'fullname': 'Peptide', 'group' : group_name, 'color' : '#666666', 'css' : css_block });
            var peptide = peps[i].sequence;
            var peptide_bits = sequenceRenderer.getAminoAcidsByPeptide(peptide);
            if (peptide_bits.length === 0){
                continue;
            }
            peptide_bits.addToLayer(layer_name);
            icons.push(peptide_bits.addToLayer(layer_name));

            for (var k = 0; k < peps[i].positions.length; k++ ) {
                icons = icons.concat(peptide_bits[peps[i].positions[k] - 1].addToLayer(overlay_name, {'content' : 'MAT', 'height' : 20, 'offset' : -2.5 }));
                peptide_bits[peps[i].positions[k] - 1].addToLayer(layer_name, {'content' : 'MAT', 'height' : 20, 'offset' : -2.5 });
            }
        }
        jQuery(sequenceRenderer).trigger('resultsRendered',[reader]);
    });
    return this;
};
/** Retrieve an array of positions that nterm has been experimentally verified to occur upon
 *  @returns {Array}    NTerm positions upon the full protein
 */
MASCP.NTermReader.Result.prototype.getAllExperimentalPositions = function()
{
    var peps = this.getPeptides();
    var results = [];
    var seen = {};
    peps.forEach(function(pep) {
        pep.positions.forEach(function(pos) {
            if ( ! seen[pos] ) {
                results.push(pos);
                seen[pos] = true;
            }
        });
    });
    return results;
}
MASCP.NTermReader.Result.prototype.render = function()
{
};
