#!/usr/bin/env python
from __future__ import with_statement

import sys
if '%x' % sys.hexversion < '2070500':
  maj, min, micro = sys.version_info[:3]
  print "Your python version '%s.%s.%s' is old, please update to the latest '2.7.x' one from 'http://www.python.org/download/'." % (maj, min,micro)
  sys.exit(1)

import os, urllib2, logging, time, urllib, socket
from argparse import ArgumentParser
from xml.dom.minidom import parse

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))
import html5lib
from bs4 import BeautifulSoup

DEFAULT_TIMEOUT = 60 * 10
PFAM_DOMAIN = 'http://pfam.sanger.ac.uk'
##PFAM_DOMAIN = 'http://pfam.janelia.org'

PUBSEED_URL = 'http://pubseed.theseed.org'


__doc__ = """This program  takes as input a txt file with a list of query protein sequences in fasta format.
For each input sequence the code identifies the corresponding Pfam protein family and queries
its "Domain organization" data.  The output file includes the list and descriptions of all fusion events
("architectures") this family is involved in. A single representative protein ID for each type of fusion events is listed.
In addition, the complete list of all fusion proteins of each architecture is provided in a separate file.
"""
program_name = os.path.basename(os.path.abspath(__file__))


def setup_logger(level):
  ### set up root level logger first
  root_logger = logging.getLogger()
  root_logger.setLevel(level)
  hdl = logging.StreamHandler()
  hdl.setFormatter(logging.Formatter('%(asctime)s %(name)s %(levelname)-8s %(message)s', '%m-%d-%y %H:%M'))
  root_logger.addHandler(hdl)


class BaseWithLogger(object):
  def __init__(self, class_name):
    self.logger = logging.getLogger(class_name)


class Member(object):
  __slots__ = ('protein_id', 'organism_name', 'function', 'length')
  def __init__(self, protein_id, organism_name, function, length):
    self.protein_id = protein_id
    self.organism_name = organism_name
    self.function = function
    self.length = length


class Architecture(BaseWithLogger):
  def __init__(self, orf_name, pfam_family, pfam_family_name, primary_id, description):
    super(Architecture, self).__init__(self.__class__.__name__)
    self.orf_name = orf_name
    self.pfam_family = pfam_family
    self.pfam_family_name = pfam_family_name
    self.primary_protein_id = primary_id
    self.description = " ".join(description.split())

    self.accessor = None
    self.members = []


  def add_member(self, text):
    protein_id = text.split('[')[0].strip()
    organism_name = text.split('[')[1].split(']')[0]
    function  = text.split(']')[1].split('(')[0].strip()
    length  = text.split('(')[-1].split('residues)')[0]

    self.members.append(Member(protein_id, organism_name, function, length))


class FastaReader(BaseWithLogger):
  def __init__(self, fname):
    super(FastaReader, self).__init__(self.__class__.__name__)

    if not os.path.exists(fname):
      self.logger.error("Could not find '%s' input fasta file.", fname)
      sys.exit(1)

    self.logger.debug("Found '%s' input fasta file.", fname)
    self.fname = fname


  @property
  def records(self):
    name, alias, seq = None, '', ''
    with open(self.fname) as f:
      current_name, current_alias = None, None
      for ln in f:
        ln = ln.strip()
        if not len(ln):
          continue

        if ln.startswith('>'):
          fields = ln.lstrip('>').split(' ')
          if len(fields) == 1:
            name = fields[0]
            alias = ''
          elif len(fields) >= 2:
            name = fields[0]
            alias = ' '.join(fields[1:])
          else:
            self.logger.error("Could not find the sequence id in '%s' fasta file", self.fname)
            sys.exit(1)

          if current_name is not None:
            yield (current_name, current_alias, seq)
            current_name = name
            current_alias = alias
            seq = ''
          else:
            current_name = name
            current_alias = alias
        else:
          seq += ln
    yield (name, alias, seq)


class PubSEEDConnector(BaseWithLogger):
  def __init__(self):
    super(PubSEEDConnector, self).__init__(self .__class__.__name__)

  #curl -d "act=do_search&page=Find&pattern=Q04CD7&submit=Search"   'http://pubseed.theseed.org' > /tmp/o
  def get_fig_id(self, aliases):
    """
    returns a tuple of the first found fig id and an alias for which it was found from the given list of aliases,
    or empty tuple, if cannot find any
    """
    for alias in aliases:
      req = urllib2.Request(PUBSEED_URL)
      req.add_header('Accept', 'text/javascript, text/html, application/xml, text/xml, */*')
      req.add_data(urllib.urlencode({'act':'do_search', 'page':'Find', 'pattern':alias.split('_')[0], 'submit':'Search'}))
      self.logger.debug("Looking for FIG ID for '%s'.", alias)
      try:
        result = urllib2.urlopen(req, timeout = DEFAULT_TIMEOUT)
      except (urllib2.URLError, socket.timeout) as e:
        if hasattr(e, 'reason'):
          reason = e.reason
        else:
          reason = e
        self.logger.error("Could not communicate with '%s', reason: '%s', please retry.", PUBSEED_URL, reason)
        sys.exit(1)

      content = result.read()
      result.close()
      soup = BeautifulSoup(content)
      for _input in soup.find_all('input'):
        if _input.has_attr('id') and _input['id']=='table_data_0':
          soup2 = BeautifulSoup(_input['value'])
          fig_id = soup2.a.text
          if len(fig_id):
            self.logger.debug("Found FIG ID '%s'." % fig_id)
            return alias,fig_id

    return '',''



class PfamConnector(BaseWithLogger):
  PFAM_FAMILY_URL = '/'.join((PFAM_DOMAIN, 'family'))
  PFAM_SEARCH_URL = '/'.join((PFAM_DOMAIN, 'search', 'sequence'))
  PFAM_DOMAIGRAPHS_URL = '/'.join((PFAM_DOMAIN,'domaingraphics'))

  def __init__(self):
    super(PfamConnector, self).__init__(self .__class__.__name__)
    self.job_ids = {}
    self.pfam_families_for_orf = {}


  def submit_sequence_for_search(self, orf_name, seq):
    req = urllib2.Request(self.PFAM_SEARCH_URL)
    req.add_data(urllib.urlencode({'seq':seq, 'output':'xml'}))
    req.add_header('Expect', '')
    self.logger.info("Submitting orf '%s' sequence for sequence search.", orf_name)
    try:
      dom = parse(urllib2.urlopen(req, timeout=DEFAULT_TIMEOUT))
    except (urllib2.URLError, socket.timeout) as e:
      if hasattr(e, 'reason'):
        reason = e.reason
      else:
        reason = e
      self.logger.error("Could not communicate with '%s', reason: '%s', please retry.", PFAM_DOMAIN, reason)
      sys.exit(1)

    try:
      ### u'http://pfam.sanger.ac.uk/search/sequence/resultset/A02BCFCC-0A0D-11E3-995B-6CD0AFE119EF?output=xml'
      job_id_url = dom.getElementsByTagName('result_url')[0].childNodes[0].data
      self.job_ids[orf_name] = job_id_url
    except IndexError:
      self.logger.error("Search server did not accept '%s' orf's sequence for search, "
                        "please try this sequence again later.", orf_name)


  def collect_results_for_search_jobs(self):
    while True:
      time.sleep(3)
      keys = self.job_ids.keys()
      if not len(keys):
        break
      for orf_name in keys:
        self.logger.info("Checking job status for orf '%s' ...", orf_name)
        jobs_id = self.job_ids.get(orf_name)
        result = urllib2.urlopen(jobs_id)
        if result.getcode() == 200:
          self.logger.debug("Received results for orf '%s'", orf_name)
          self.__get_pfam_families_from_search_results(orf_name, result)
          self.job_ids.pop(orf_name)
        elif result.getcode() == 202:
          self.logger.info("Still waiting on results for orf '%s'", orf_name)
        elif result.getcode() == 502:
          self.logger.info("The job for orf '%s' failed on the search system side. "
                           "Please try to resubmit this orf again.", orf_name)
          self.job_ids.pop(orf_name)
        elif result.getcode() == 503:
          self.logger.info("The job for orf '%s' is put on hold by server Admin. "
                           "You may want to contact them or wait and try to resubmit this orf again.", orf_name)
          self.job_ids.pop(orf_name)
        elif result.getcode() == 401:
          self.logger.info("The job for orf '%s' was deleted from the search system by server Admin. "
                           "There was probably a problem with the job and you should contact the help "
                           "desk for assistance with it, or wait and try to resubmit this orf again.", orf_name)
          self.job_ids.pop(orf_name)
        elif result.getcode() == 500:
          self.logger.info("The job for orf '%s' is put on hold by server Admin. "
                           "You may want to contact them or wait and try to resubmit this orf again.", orf_name)
          self.job_ids.pop(orf_name)
        else:
          self.logger.warning("Got unexpected status code '%s' for orf '%s'. Do not know what to do, will skip "
                              "this orf. Please contact developer.", result.getcode(), orf_name)
          self.job_ids.pop(orf_name)


  def __get_pfam_families_from_search_results(self, orf_name, result_object):
    dom = parse(result_object)
    for element in dom.getElementsByTagName('match'):
      pfam_family = element.attributes.getNamedItem('accession').value
      self.logger.debug("For orf '%s' found pfam family '%s'.", orf_name, pfam_family)
      self.pfam_families_for_orf.setdefault(orf_name,[]).append(pfam_family)


  def __get_pfam_family_name(self, pfam_family):
    url = '/'.join((self.PFAM_FAMILY_URL,pfam_family)) + '?output=xml'
    self.logger.debug("Retrieving family name for '%s'.", pfam_family)
    soup = BeautifulSoup(self.__get_content(url))
    return soup.entry['id']


  def __get_content(self, url):
    req = urllib2.Request(url)
    req.add_header('Expect', '')
    req.add_header('Accept', 'text/javascript, text/html, application/xml, text/xml, */*')
    ##req.add_header('Accept-Encoding',	'gzip, deflate')
    ##req.add_header('User-Agent', 'curl/7.30.0')
    ##req.add_header('User-Agent', 'curl/7.24.0 (x86_64-apple-darwin12.0) libcurl/7.24.0 OpenSSL/0.9.8x zlib/1.2.5')
    self.logger.debug("Retrieving content from '%s'...", url)
    try:
      result = urllib2.urlopen(req, timeout = DEFAULT_TIMEOUT)
    except (urllib2.URLError, socket.timeout) as e:
      if hasattr(e, 'reason'):
        reason = e.reason
      else:
        reason = e
      self.logger.error("Could not communicate with '%s', reason: '%s', please retry.", PFAM_DOMAIN, reason)
      sys.exit(1)

    content = result.read()
    result.close()
    return content


  def process_found_pfam_families(self):
    all_architectures = []
    for orf_name, pfam_families in self.pfam_families_for_orf.iteritems():
      self.logger.info("Processing pfam families for orf '%s' ...", orf_name)
      for pfam_family in pfam_families:
        architectures = []
        self.logger.debug("Processing pfam family '%s'.", pfam_family)

        pfam_family_name = self.__get_pfam_family_name(pfam_family)

        self.logger.debug("Requesting architectures for pfam family '%s'.", pfam_family)
        url = '/'.join((self.PFAM_DOMAIGRAPHS_URL, pfam_family))
        soup = BeautifulSoup(self.__get_content(url), 'html5lib')

        for div in soup.find_all('div'):
          if div.has_attr('class') and div.has_attr('id') and ('graphicRow' in div['class']) and div['id'].startswith('row'):
              primary_id = div['id'].lstrip('row')
              description = div.h3.text.strip()
              architecture_name = description.split('architecture:')[1]

              ## keep lines with ' x ' or ','
              if (' x ' in architecture_name) or (',' in architecture_name):
                self.logger.debug("Adding architecture with primary protein id '%s'.", primary_id)
                architectures.append(Architecture(orf_name,pfam_family,pfam_family_name,primary_id,description))

        lines = str(soup.script.string.split('var layout')[0]).split('\n')
        for n, ln in enumerate(lines):
          ln = ln.strip()
          for architecture in architectures:
            if 'row%s' % architecture.primary_protein_id in ln:
              accessor = lines[n+1].split('.store( "arch",')[1].strip().replace('"','').replace(' );','')
              architecture.accessor = accessor

              self.logger.debug("Requesting architecture details for accessor '%s', may take some time, please wait ...", accessor)
              url = '/'.join((self.PFAM_DOMAIGRAPHS_URL,architecture.pfam_family))
              url = '?arch='.join((url,accessor))
              soup = BeautifulSoup(self.__get_content(url))

              for item in soup.div.find_all('div'):
                if not len(item):
                  continue
                text = ' '.join(item.text.split())
                architecture.add_member(text)

        all_architectures.extend(architectures)
    return all_architectures


def process(options):
  setup_logger(options.vlevel)
  logger = logging.getLogger(program_name)

  orf_name_to_alias = {}
  fasta_reader = FastaReader(options.input_file[0])
  pfam_connector = PfamConnector()
  for orf_name, alias, seq in fasta_reader.records:
    logger.debug("Read fasta record for ORF '%s'.", orf_name)
    orf_name_to_alias[orf_name] = alias
    pfam_connector.submit_sequence_for_search(orf_name, seq)
    time.sleep(2)

  pfam_connector.collect_results_for_search_jobs()
  architectures = pfam_connector.process_found_pfam_families()

  seed_connector = PubSEEDConnector()
  logger.info("Looking for FIG ID's ...")
  with open(options.output_summary_file[0],'w') as f:
    header = "ORF ID\tFunction\tPfam family\tNumber of seqs with this architecture\tRepresentative Protein ID\tFIG ID\tDomains"
    f.write("%s\n" % header)

    for architecture in architectures:
      if 'There is' in architecture.description:
        num_of_seq = architecture.description.split('There is')[1].strip().split()[0]
      else:
        num_of_seq = architecture.description.split('There are')[1].strip().split()[0]
      domains = '\t'.join(architecture.description.split('architecture:')[1].strip().split(','))

      protein_id, fig_id = seed_connector.get_fig_id([member.protein_id for member in architecture.members])
      if len(fig_id):
        fig_id = '=HYPERLINK("%s/?page=Annotation&feature=%s","%s")' % (PUBSEED_URL, fig_id, fig_id)
      else:
        protein_id = architecture.primary_protein_id

      f.write('%s\t%s\t%s (%s)\t%s\t=HYPERLINK("%s/protein/%s","%s")\t%s\t%s\n'
             %(architecture.orf_name, orf_name_to_alias[architecture.orf_name],
             architecture.pfam_family_name, architecture.pfam_family,
             num_of_seq, PFAM_DOMAIN, protein_id, protein_id, fig_id, domains))


  with open(options.output_details_file[0],'w') as f:
    header = 'ORF ID\tFunction\tPfam family\tArchitecture, with statistics\tProtein ID\tFunction\tLength\tOrganism Name'
    f.write("%s\n" % header)

    for architecture in architectures:
      for member in architecture.members:
        f.write('%s\t%s\t%s (%s)\t%s\t=HYPERLINK("%s/protein/%s","%s")\t%s\t%s\t%s\n'
              %(architecture.orf_name, orf_name_to_alias[architecture.orf_name],
                architecture.pfam_family_name, architecture.pfam_family, architecture.description,
                PFAM_DOMAIN, member.protein_id, member.protein_id,
                member.function, member.length, member.organism_name))

  return 0


def main():
  parser = ArgumentParser(description = __doc__)
  parser.add_argument('-v', '--vlevel', type=int, default=20, help="Verbosity logging level. "
                      "Available levels are:  50:CRITICAL; 40:ERROR; 30:WARNING; 20:INFO; 10:DEBUG. Default: 20")
  parser.add_argument('input_file', nargs=1, help='Fasta formatted input file name.' )
  parser.add_argument('output_summary_file', nargs=1, help='The summary result file name.')
  parser.add_argument('output_details_file', nargs=1, help='The details result file name.')
  options = parser.parse_args()
  return process(options)

if __name__ == "__main__":
  sys.exit(main())
