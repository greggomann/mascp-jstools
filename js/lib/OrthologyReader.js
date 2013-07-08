/** @fileOverview   Classes for reading data from the PPAPR orthology database
 */
if ( typeof MASCP == 'undefined' || typeof MASCP.Service == 'undefined' ) {
    throw "MASCP.Service is not defined, required class";
}

/** Default class constructor
 *  @class      Service class that will retrieve data from Orthology for a given AGI.
 *  @param      {String} agi            Agi to look up
 *  @param      {String} endpointURL    Endpoint URL for this service
 *  @extends    MASCP.Service
 */
MASCP.OrthologyReader = MASCP.buildService(function(data) {
                        this._raw_data = data;                        
                        return this;
                    });

MASCP.OrthologyReader.SERVICE_URL = 'http://prabi2.inrialpes.fr/at_chloro/annotation/';

MASCP.OrthologyReader.prototype.requestData = function()
{
    var agi = this.agi;
    
    return {
        type: "GET",
        url: this._endpointURL + agi.toUpperCase(),
        dataType: "json",
        data: { 'agi'       : agi.toUpperCase(),
                'service'   : 'orthology' 
        }
    };
};


/**
 *  @class   Container class for results from the Orthology service
 *  @extends MASCP.Service.Result
 */
// We need this line for the JsDoc to pick up this class
MASCP.OrthologyReader.Result = MASCP.OrthologyReader.Result;

/** Retrieve the peptides for this particular entry from the Orthology service
 *  @returns Array of peptide strings
 *  @type [String]
 */
MASCP.OrthologyReader.Result.prototype.getConservation = function()
{
    var content = null;

    if (this._conservation) {
        return this._conservation;
    }

    this._long_name_map = {};
    
    if (! this._raw_data || ! this._raw_data.conservation ) {
        return [];
    }

        
    var conservation = [];
    
    for (var i = 0; i < this._raw_data.conservation.length; i++ ) {
        var a_conservation = this._raw_data.conservation[i];
        conservation.push(a_conservation);
    }
    this._conservation = conservation;
    return conservation;
};

MASCP.OrthologyReader.Result.prototype._cleanSequence = function(sequence)
{
    return sequence.replace(/[^A-Z]/g,'');
};

MASCP.OrthologyReader.prototype.setupSequenceRenderer = function(sequenceRenderer)
{
    var reader = this;

    var css_block = '.active .overlay { background: #55ff33; } .active a { color: #000000; text-decoration: none !important; }  :indeterminate { background: #ff0000; } .tracks .active { background: #0000ff; } .inactive a { text-decoration: none; } .inactive { display: none; }';
    

    this.bind('resultReceived', function() {
        var conservation = this.result.getConservation();
        if (conservation.length > 0) {
            jQuery(this.renderers).each(function(i){
                this.createOrthology(conservation);
            });
            MASCP.getLayer('orthology_experimental').href = 'http://prabi2.inrialpes.fr/at_chloro/protein/'+reader.agi.toUpperCase();
        }
        jQuery(sequenceRenderer).trigger('resultsRendered',[reader]);
    });
    return this;
};

MASCP.OrthologyReader.Result.prototype.render = function()
{
    if (this.getConservation().length > 0) {
        var a_container = jQuery('<div>MS/MS spectra <input class="group_toggle" type="checkbox"/>Orthology</div>');
        jQuery(this.reader.renderers).each(function(i){
            this.createGroupCheckbox('orthology_experimental',jQuery('input.group_toggle',a_container));
        });
        return a_container;
    } else {
        return null;
    }
};
