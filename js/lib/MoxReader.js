/** @fileOverview   Classes for reading data from the Mox data
 */
if ( typeof MASCP == 'undefined' || typeof MASCP.Service == 'undefined' ) {
    throw "MASCP.Service is not defined, required class";
}

/** Default class constructor
 *  @class      Service class that will retrieve data from the Mox data for a given AGI.
 *  @param      {String} agi            Agi to look up
 *  @param      {String} endpointURL    Endpoint URL for this service
 *  @extends    MASCP.Service
 */
MASCP.MoxReader = MASCP.buildService(function(data) {
                        this._raw_data = data;
                        return this;
                    });

MASCP.MoxReader.SERVICE_URL = '?';

MASCP.MoxReader.prototype.requestData = function()
{
    var agi = this.agi;

    return {
        type: "GET",
        dataType: "json",
        data: { 'agi'       : agi,
                'service'   : 'mox'
        }
    };
};

/**
 *  @class   Container class for results from the Mox service
 *  @extends MASCP.Service.Result
 */
// We need this line for the JsDoc to pick up this class
MASCP.MoxReader.Result = MASCP.MoxReader.Result;

/** Retrieve the peptides for this particular entry from the Mox service
 *  @returns Array of peptide strings
 *  @type [String]
 */
MASCP.MoxReader.Result.prototype.getPeptides = function()
{
    var content = null;
    if (! this._raw_data || ! this._raw_data.data  || ! this._raw_data.data.peptides ) {
        return [];
    }

    return this._raw_data.data.peptides;
};

MASCP.MoxReader.Result.prototype._cleanSequence = function(sequence)
{
    return sequence.replace(/[^A-Z]/g,'');
};

MASCP.MoxReader.prototype.setupSequenceRenderer = function(sequenceRenderer)
{
    var reader = this;

    var css_block = '.active .overlay { background: #666666; } .active a { color: #000000; text-decoration: none !important; }  :indeterminate { background: #ff0000; } .tracks .active { background: #0000ff; } .inactive a { text-decoration: none; } .inactive { display: none; }';

    this.bind('resultReceived', function() {
        var peps = this.result.getPeptides();

        var overlay_name = 'mox_experimental';
        var group_name = 'mox_peptides';
        var icons = [];

        if (peps.length > 0) {
            MASCP.registerLayer(overlay_name,{ 'fullname' : 'Met Ox (mod)', 'color' : '#666666', 'css' : css_block });

            MASCP.registerGroup(group_name, {'fullname' : 'Met Ox (mod)', 'hide_member_controllers' : true, 'hide_group_controller' : true, 'color' : '#666666' });
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
            var layer_name = 'mox_peptide_'+i;
            MASCP.registerLayer(layer_name, { 'fullname': 'Peptide', 'group' : group_name, 'color' : '#666666', 'css' : css_block });
            var peptide = peps[i].sequence;
            var peptide_bits = sequenceRenderer.getAminoAcidsByPeptide(peptide);
            if (peptide_bits.length === 0){
                continue;
            }
            peptide_bits.addToLayer(layer_name);
            icons.push(peptide_bits.addToLayer(layer_name));

            for (var k = 0; k < peps[i].positions.length; k++ ) {
                icons = icons.concat(peptide_bits[peps[i].positions[k] - 1].addToLayer(overlay_name, {'content' : 'MOX', 'height' : 20, 'offset' : -2.5 }));
                peptide_bits[peps[i].positions[k] - 1].addToLayer(layer_name, {'content' : 'MOX', 'height' : 20, 'offset' : -2.5 });
            }
        }
        jQuery(sequenceRenderer).trigger('resultsRendered',[reader]);
    });
    return this;
};
/** Retrieve an array of positions that mox has been experimentally verified to occur upon
 *  @returns {Array}    Mox positions upon the full protein
 */
MASCP.MoxReader.Result.prototype.getAllExperimentalPositions = function()
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
MASCP.MoxReader.Result.prototype.render = function()
{
};
