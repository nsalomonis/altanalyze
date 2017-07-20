import base64
import csv
import getpass
try:
    import lxml
except Exception:
    print 'library lxml not supported. WikiPathways and LineageProfiler visualization will not work. Please install with pip install lxml.'
from lxml import etree as ET
from lxml import _elementpath
import re
try: import requests
except Exception:
    print 'library requests not supported. WikiPathways and LineageProfiler visualization will not work. Please install with pip install requests.'
import sys

class WikipathwaysApiClient(object):
    """Returns :class:`WikipathwaysApiClient` object.
    :param identifier: WikiPathways ID for the new :class:`WikipathwaysApiClient` object.
    """

    def __invert_dict(dictionary):
        return dict((v, k) for k, v in dictionary.iteritems())

    def __get_bridgedb_datasets(self):
        if hasattr(self, 'bridgedb_datasets'):
            bridgedb_datasets = self.bridgedb_datasets
        else:
            bridgedb_datasets_request = requests.get('https://raw.githubusercontent.com/bridgedb/BridgeDb/master/org.bridgedb.bio/resources/org/bridgedb/bio/datasources.txt')
            bridgedb_datasets_string = bridgedb_datasets_request.text
            bridgedb_datasets_csv = csv.reader(bridgedb_datasets_string.split('\n'), delimiter='\t')
            
            bridgedb_datasets_parsed = [];
            for bridgedb_dataset_csv in bridgedb_datasets_csv:
                if bridgedb_dataset_csv:
                    bridgedb_dataset_parsed = {}
                    bridgedb_dataset_parsed['system_code'] = bridgedb_dataset_csv[1]
                    bridgedb_dataset_parsed['miriam'] = bridgedb_dataset_csv[8]
                    bridgedb_datasets_parsed.append(bridgedb_dataset_parsed)

            self.bridgedb_datasets = bridgedb_datasets_parsed

        return self.bridgedb_datasets

    def __parse_identifiers_iri(self, iri):
        iri_components = iri.split('identifiers.org')
        iri_path = iri_components[len(iri_components) - 1]
        iri_path_components = iri_path.split('/')
        preferred_prefix = iri_path_components[1]
        identifier = iri_path_components[2]

        bridgedb_datasets = self.__get_bridgedb_datasets()
        for bridgedb_dataset in bridgedb_datasets:
            if 'urn:miriam:' + preferred_prefix == bridgedb_dataset['miriam']:
                system_code = bridgedb_dataset['system_code']

        return {
            'preferred_prefix': preferred_prefix,
            'identifier': identifier,
            'system_code': system_code
        }

    api_to_standard_term_mappings = {
        'id': 'identifier',
        'ids': 'identifiers',
        'pwId': 'identifier',
        'revision': 'version',
        'graphId': 'element_identifiers',
        'color': 'colors',
        'fileType': 'file_format',
        'species': 'organism',
        'url': 'web_page',
        'codes': 'system_codes'
    }

    
    filename_extension_to_media_type_mappings = {
        'svg': 'image/svg+xml',
        'png': 'image/png',
        'pdf': 'application/pdf',
        'gpml': 'application/gpml+xml',
        'txt': 'text/vnd.genelist+tab-separated-values',
        'pwf': 'text/vnd.eu.gene+plain',
        'owl': 'application/vnd.biopax.owl+xml',
    }
    filename_extensions = filename_extension_to_media_type_mappings.keys()
    media_types = filename_extension_to_media_type_mappings.values()
    media_type_to_filename_extension_mappings = __invert_dict(filename_extension_to_media_type_mappings)


    english_name_to_iri_mappings = {
            'African malaria mosquito': 'http://identifiers.org/taxonomy/7165',
            'beet': 'http://identifiers.org/taxonomy/161934',
            'thale cress': 'http://identifiers.org/taxonomy/3702',
            'cattle': 'http://identifiers.org/taxonomy/9913',
            'roundworm': 'http://identifiers.org/taxonomy/6239',
            'dog': 'http://identifiers.org/taxonomy/9615',
            'sea vase': 'http://identifiers.org/taxonomy/7719',
            'zebrafish': 'http://identifiers.org/taxonomy/7955',
            'fruit fly': 'http://identifiers.org/taxonomy/7227',
            'Escherichia coli': 'http://identifiers.org/taxonomy/562',
            'horse': 'http://identifiers.org/taxonomy/9796',
            'chicken': 'http://identifiers.org/taxonomy/9031',
            'soybean': 'http://identifiers.org/taxonomy/3847',
            'human': 'http://identifiers.org/taxonomy/9606',
            'barley': 'http://identifiers.org/taxonomy/4513',
            'Rhesus monkey': 'http://identifiers.org/taxonomy/9544',
            'mouse': 'http://identifiers.org/taxonomy/10090',
            'platypus': 'http://identifiers.org/taxonomy/9258',
            'long-grained rice': 'http://identifiers.org/taxonomy/39946',
            'rice': 'http://identifiers.org/taxonomy/4530',
            'black cottonwood': 'http://identifiers.org/taxonomy/3694',
            'chimpanzee': 'http://identifiers.org/taxonomy/9598',
            'Norway rat': 'http://identifiers.org/taxonomy/10116',
            'baker\'s yeast': 'http://identifiers.org/taxonomy/4932',
            'tomato': 'http://identifiers.org/taxonomy/4081',
            'pig': 'http://identifiers.org/taxonomy/9823',
            'wine grape': 'http://identifiers.org/taxonomy/29760',
            'western clawed frog': 'http://identifiers.org/taxonomy/8364',
            'maize': 'http://identifiers.org/taxonomy/4577'
    }
    iri_to_english_name_mappings = __invert_dict(english_name_to_iri_mappings)


    latin_name_to_iri_mappings = {
            'Anopheles gambiae': 'http://identifiers.org/taxonomy/7165',
            'Arabidopsis thaliana': 'http://identifiers.org/taxonomy/3702',
            'Aspergillus niger': 'http://identifiers.org/taxonomy/5061',
            'Bacillus subtilis': 'http://identifiers.org/taxonomy/1423',
            'Beta vulgaris': 'http://identifiers.org/taxonomy/161934',
            'Bos taurus': 'http://identifiers.org/taxonomy/9913',
            'Caenorhabditis elegans': 'http://identifiers.org/taxonomy/6239',
            'Canis familiaris': 'http://identifiers.org/taxonomy/9615',
            'Ciona intestinalis': 'http://identifiers.org/taxonomy/7719',
            'Ruminiclostridium thermocellum': 'http://identifiers.org/taxonomy/1515',
            'Clostridium thermocellum': 'http://identifiers.org/taxonomy/1515',
            'Danio rerio': 'http://identifiers.org/taxonomy/7955',
            'Drosophila melanogaster': 'http://identifiers.org/taxonomy/7227',
            'Escherichia coli': 'http://identifiers.org/taxonomy/562',
            'Equus caballus': 'http://identifiers.org/taxonomy/9796',
            'Gallus gallus': 'http://identifiers.org/taxonomy/9031',
            'Gibberella zeae': 'http://identifiers.org/taxonomy/5518',
            'Glycine max': 'http://identifiers.org/taxonomy/3847',
            'Homo sapiens': 'http://identifiers.org/taxonomy/9606',
            'Hordeum vulgare': 'http://identifiers.org/taxonomy/4513',
            'Macaca mulatta': 'http://identifiers.org/taxonomy/9544',
            'Mus musculus': 'http://identifiers.org/taxonomy/10090',
            'Mycobacterium tuberculosis': 'http://identifiers.org/taxonomy/1773',
            'Ornithorhynchus anatinus': 'http://identifiers.org/taxonomy/9258',
            'Oryza indica': 'http://identifiers.org/taxonomy/39946',
            'Oryza sativa': 'http://identifiers.org/taxonomy/4530',
            'Oryza sativa Indica Group': 'http://identifiers.org/taxonomy/39946',
            'Populus trichocarpa': 'http://identifiers.org/taxonomy/3694',
            'Pan troglodytes': 'http://identifiers.org/taxonomy/9598',
            'Rattus norvegicus': 'http://identifiers.org/taxonomy/10116',
            'Saccharomyces cerevisiae': 'http://identifiers.org/taxonomy/4932',
            'Solanum lycopersicum': 'http://identifiers.org/taxonomy/4081',
            'Sus scrofa': 'http://identifiers.org/taxonomy/9823',
            'Vitis vinifera': 'http://identifiers.org/taxonomy/29760',
            'Xenopus tropicalis': 'http://identifiers.org/taxonomy/8364',
            'Zea mays': 'http://identifiers.org/taxonomy/4577'
    }


    def __init__(self, base_iri=None):
        if base_iri is None:
            base_iri = 'http://webservice.wikipathways.org/'
        self.base_iri = base_iri

        # define namespaces
        self.NAMESPACES = {'ns1':'http://www.wso2.org/php/xsd','ns2':'http://www.wikipathways.org/webservice/'}


    def __convert_standard_terms_to_api_terms(self, input_params):
        terms_to_convert = self.api_to_standard_term_mappings
        standard_terms = terms_to_convert.values()
        api_terms = terms_to_convert.keys()

        request_params = {}
        for key, value in input_params.iteritems():
            if (key in standard_terms):
                def get_api_term(candidate_api_term):
                    return self.api_to_standard_term_mappings[candidate_api_term] == key

                api_term = filter(get_api_term, api_terms)[0]
                request_params[api_term] = input_params[key]
            else:
                request_params[key] = input_params[key]

        return request_params


    def __convert_api_terms_to_standard_terms(self, input_object):
        terms_to_convert = self.api_to_standard_term_mappings
        standard_terms = terms_to_convert.values()
        api_terms = terms_to_convert.keys()

        output_object = {}
        for key, value in input_object.iteritems():
            if (key in api_terms):
                api_term = terms_to_convert[key]
                output_object[api_term] = input_object[key]
            else:
                output_object[key] = input_object[key]

        return output_object


    def __convert_organism_to_dict(self, organism):
        if hasattr(self, 'organism_dicts'):
            organism_dicts = self.organism_dicts
        else:
            organism_dicts = self.organism_dicts = self.list_organisms()

        for organism_dict in organism_dicts:
            if isinstance(organism, basestring):
                if organism_dict['@id'] == organism:
                    return organism_dict
                elif organism_dict['name']['la'] == organism:
                    return organism_dict
                elif organism_dict['name'].get('en') and organism_dict['name']['en'] == organism:
                    return organism_dict
            elif organism.get('@id') and organism['@id'] == organism_dict['@id']:
                    return organism_dict


    def __enrich_pathway(self, pathway):
        pathway['@id'] = 'http://identifiers.org/wikipathways/' + pathway['identifier']

        if pathway.get('organism') and isinstance(pathway['organism'], basestring):
            pathway['organism'] = self.__convert_organism_to_dict(pathway['organism'])

        return pathway


    def create_pathway(self):
        ###
        # author: msk (mkutmon@gmail.com)
        ###

        # login
        pswd = getpass.getpass('Password:')
        auth = {'name' : username , 'pass' : pswd}
        r_login = requests.get(self.base_iri + 'login', params=auth)
        dom = ET.fromstring(r_login.text)

        authentication = ''
        for node in dom.findall('ns1:auth', namespaces):
                authentication = node.text

        # read gpml file
        f = open(gpml_file, 'r')
        gpml = f.read()

        # create pathway
        update_params = {'auth' : username+'-'+authentication, 'gpml': gpml}
        re = requests.post(self.base_iri + 'createPathway', params=update_params)
        #print re.text


    def get_colored_pathway(self, identifier, element_identifiers, colors, version = '0', file_format = 'svg'):
        """Sends a GET request. Returns file as string.

        Args:
          identifier (str): WikiPathways ID.
          element_identifiers (list of str): means of identifying one or more elements in a pathway,
            for example, specify GPML GraphIds as ["ffffff90","ffffffe5"].
          colors (list of str): one or more hexadecimal number(s), representing the colors to use for
            the corresponding element_identifier (if the length of the colors list is equal to the
            length of the element_identifiers list) or the single color to use for all element_identifiers
            (if the colors list is not equal in length to the element_identifiers list).
            Example: ["#0000FF","#0000FF"].
          version (str, optional): The version of the pathway. Defaults to '0', which means latest.
          file_format (str): IANA media type (http://www.iana.org/assignments/media-types/media-types.xhtml)
            or filename extension desired for response. Defaults to 'svg'. Examples:

            Media types:
            * 'image/svg+xml'
            * 'image/png'
            * 'application/pdf'
            Filename extensions:
            * 'svg'
            * 'png'
            * 'pdf'
        """

        # API does not yet support content-type negotiation, so we need to convert
        # filename extension to be used as a query parameter.
        if file_format in self.media_types:
            file_format = self.media_type_to_filename_extension_mappings[file_format]

        # HTML/CSS defaults use a pound sign before the HEX number, e.g. #FFFFFF.
        # But the API does not use this, so to make it easier for users, we are
        # accepting the pound sign in the input args and stripping it here.
        input_colors = colors
        colors = []
        non_letter_number_pattern = re.compile('[^a-zA-Z0-9]+')
        for input_color in input_colors:
            color = non_letter_number_pattern.sub('', input_color)
            colors.append(color)

        input_params = {
            'identifier': identifier,
            'version': version,
            'element_identifiers': element_identifiers,
            'colors': colors,
            'file_format': file_format
        }
        request_params = self.__convert_standard_terms_to_api_terms(input_params)
        response = requests.get(self.base_iri + 'getColoredPathway', params=request_params)
        dom = ET.fromstring(response.text)
        node = dom.find('ns1:data', self.NAMESPACES)
        file = base64.b64decode(node.text) ### decode this file
        return file


    def get_pathway_as(self, identifier, version = '0', file_format = 'gpml'):
        """
        Sends a GET request. Returns an LXML object for any XML media type
        and a string for anything else

        Args:
          identifier (str): WikiPathways ID.
          version (str, optional): The version of the pathway. Defaults to '0', which means latest.
          file_format (str): IANA media type (http://www.iana.org/assignments/media-types/media-types.xhtml)
            or filename extension desired for response. Defaults to 'gpml'.
            Examples:

            Media types:
            * 'application/gpml+xml'
            * 'text/vnd.genelist+tab-separated-values'
            * 'text/vnd.eu.gene+plain'
            * 'application/vnd.biopax.owl+xml'
            * 'image/svg+xml'
            * 'image/png'
            * 'application/pdf'
            Filename extensions:
            * 'gpml'
            * 'txt'
            * 'pwf'
            * 'owl'
            * 'svg'
            * 'png'
            * 'pdf'
        """

        # API does not yet support content-type negotiation, so we need to convert
        # filename extension to be used as a query parameter.
        if file_format in self.media_types:
            file_format = self.media_type_to_filename_extension_mappings[file_format]

        input_params = {
            'identifier': identifier,
            'version': version,
            'file_format': file_format
        }

        request_params = self.__convert_standard_terms_to_api_terms(input_params)
        response = requests.get(self.base_iri + 'getPathwayAs', params=request_params)
        #print [response.text];sys.exit()
        dom = ET.fromstring(response.text)
        node = dom.find('ns1:data', self.NAMESPACES)
        response_string = base64.b64decode(node.text) ### decode this file
        if request_params['fileType'] == 'gpml' or request_params['fileType'] == 'owl' or request_params['fileType'] == 'svg':
            #response = ET.fromstring(response_string)
            response = response_string
        else:
            response = response_string
        return response


    def get_pathway_info(self, identifier):
        """Sends a GET request. Returns pathway metadata as dict.
        Args:
          identifier (str): WikiPathways ID.
        """

        request_params = {'pwId' : identifier}
        response = requests.get(self.base_iri + 'getPathwayInfo', params=request_params)
        dom = ET.fromstring(response.text)

        pathway_using_api_terms = {}
        for node in dom.findall('ns1:pathwayInfo', self.NAMESPACES):
            for attribute in node:
                pathway_using_api_terms[ET.QName(attribute).localname] = attribute.text
        pathway = self.__convert_api_terms_to_standard_terms(pathway_using_api_terms)
        pathway = self.__enrich_pathway(pathway)

        return pathway


    def find_pathways_by_text(self, query, organism = None):
        """Sends a GET request. Returns pathways as list of dicts.
        Args:
          query (str): Text to search for.
          organism (str or dict, optional): Limit to organism with given name
            (Latin or English) or @id (from http://identifiers.org/taxonomy/)
        """

        input_params = {}
        input_params['query'] = query

        if organism:
            input_params['organism'] = self.__convert_organism_to_dict(organism)['name']['la']

        request_params = self.__convert_standard_terms_to_api_terms(input_params)
        response = requests.get(self.base_iri + 'findPathwaysByText', params=request_params)
        dom = ET.fromstring(response.text)

        pathways = []
        for node in dom.findall('ns1:result', self.NAMESPACES):
            pathway_using_api_terms = {}
            for child in node:
                pathway_using_api_terms[ET.QName(child).localname] = child.text
            pathway = self.__convert_api_terms_to_standard_terms(pathway_using_api_terms)
            pathway = self.__enrich_pathway(pathway)
            pathways.append(pathway)
        return pathways


    def find_pathways_by_xref(self, **kwargs):
        """Sends a GET request. Returns pathways as a list of dicts.
            Required: either just @id or both system_codes and identifiers.
        Args:
          @id (list of str): One or more identifiers.org IRIs, like ['http://identifiers.org/ncbigene/7040'].
          system_codes (list of str): One or more BridgeDb system codes.
          identifiers (list of str): One or more entity reference identifiers.
        """

        if kwargs.get('@id'):
            if not isinstance(kwargs['@id'], list):
                kwargs['@id'] = [kwargs['@id']]
            system_codes = []
            identifiers = []
            for iri in kwargs['@id']:
                identifiers_iri_components = self.__parse_identifiers_iri(iri)
                system_codes.append(identifiers_iri_components['system_code'])
                identifiers.append(identifiers_iri_components['identifier'])
            input_params = {
                'system_codes': system_codes,
                'identifiers': identifiers
            }
        else:
            input_params = kwargs

        request_params = self.__convert_standard_terms_to_api_terms(input_params)
        response = requests.get(self.base_iri + 'findPathwaysByXref', params=request_params)
        dom = ET.fromstring(response.text)

        pathways = []
        for resultNode in dom.findall('ns1:result', self.NAMESPACES):
            pathway_using_api_terms = {}
            pathway_using_api_terms['fields'] = []
            for childNode in resultNode:
                if ET.QName(childNode).localname != 'fields':
                    pathway_using_api_terms[ET.QName(childNode).localname] = childNode.text
                elif ET.QName(childNode).localname == 'fields':
                    field = {}
                    for fieldChildNode in childNode:
                        #TODO standardize names & values from fieldChildNode.text
                        field[ET.QName(fieldChildNode).localname] = fieldChildNode.text
                    pathway_using_api_terms['fields'].append(field)
            pathway = self.__convert_api_terms_to_standard_terms(pathway_using_api_terms)
            pathway = self.__enrich_pathway(pathway)
            pathways.append(pathway)
        return pathways


    def list_organisms(self):
        """Sends a GET request. Returns :list:`organisms` object, each an organism as a dict,
        with the IRI, Latin name and English name (when available).
        """
        response = requests.get(self.base_iri + 'listOrganisms')
        dom = ET.fromstring(response.text)

        organisms = []
        for node in dom:
            try:
                organism = {}
                organism['@context'] = [
                  {
                    'name': {
                      '@id': 'biopax:name',
                      '@container': '@language'
                    },
                    'Organism': 'http://identifiers.org/snomedct/410607006'
                  }
                ]
                organism['@type'] = 'Organism'
                organism['name'] = {}
                organism['name']['la'] = latin_name = node.text
                organism['@id'] = self.latin_name_to_iri_mappings[latin_name]
                english_name = self.iri_to_english_name_mappings.get(organism['@id'])
                if english_name != None:
                    organism['name']['en'] = english_name
                organisms.append(organism)  
            except Exception:
                pass


        return organisms


    # list pathways
    def list_pathways(self, organism):
        request_params = {'organism': organism}
        response = requests.get(self.base_iri + 'listPathways', params=request_params)
        dom = ET.fromstring(response.text)

        pathways = []
        for pathway_node in dom.findall('ns1:pathways', self.NAMESPACES):
            pathway_using_api_terms = {}
            for child_node in pathway_node:
                pathway_using_api_terms[ET.QName(child_node).localname] = child_node.text
            pathway = self.__convert_api_terms_to_standard_terms(pathway_using_api_terms)
            pathway = self.__enrich_pathway(pathway)
            pathways.append(pathway)

        return pathways

if __name__ == '__main__':
    client = WikipathwaysApiClient()
    wp_id_data = client.get_pathway_as(file_format = 'gpml',identifier = 'WP254', version = 0)
    with open('WP205.gpml', 'a') as file_out:
        file_out.write(wp_id_data)