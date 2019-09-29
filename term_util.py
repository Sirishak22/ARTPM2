##############################################
# Term manipulation and Indri query language methods.
# Lowell Milliken
##############################################
import re


def clean(s):
    """Remove non-alpanumerics."""
    notalphanum = r'[^a-zA-Z0-9_ ]'
    return re.sub(notalphanum, ' ', s)


def exact(s):
    """Add exact phrase Indri tag."""
    if s:
        return '#1({})'.format(s)
    else:
        return ''


def combine(s):
    """Add combine Indri tage"""
    return '#combine({})'.format(s)


def field(s, field):
    """Add term field restrictions"""
    terms = s.split()
    field_s = ''
    for term in terms:
        field_s += ' ' + '{}.{}'.format(term, field)
    return field_s


def combine_field(s, field):
    """Add combine field restrictions"""
    return '#combine[{}]({})'.format(field, s)


def form_query(number, s, pmids=None):
    """Put the query string in XML tags."""
    if pmids is not None:
        query = '<query><number>{}</number><text>{}</text>'.format(number, s)
        query += ''.join(['<workingSetDocno>{}</workingSetDocno>'.format(pmid) for pmid in pmids[number]])
        query += '</query>\n'
        return query
    else:
        return '<query><number>{}</number><text>{}</text></query>\n'.format(number, s)


def syn(s):
    """Add synonyms tag."""
    if s.strip():
        return '#syn({})'.format(s)
    return s


def band(s1, s2):
    """Add band between two strings."""
    return '#band({} {})'.format(s1, s2)


def scoreif(f, s):
    """Add scoreif to the beginning of a query string.

    :param f: filter string
    :param s: query string
    :return:
    """
    return '#scoreif({} {})'.format(f, s)


def uwindow(s):
    """Put terms in an unordered window."""
    return '#uw({})'.format(s)


def find_age_group(demographic):
    """Get age group words out of the raw demographic."""
    age_groups = [['pediatric', 'paediatric', 'children'],
                  ['adult'],
                  ['elderly', 'old', 'geriatric', 'aged', 'senescence', 'senium']]
    age_thresholds = [19, 65]
    age = int(demographic.split(' ')[0])
    for i in range(len(age_thresholds)):
        if age < age_thresholds[i]:
            return age_groups[i]  # , age_group_cuis[i]

    return age_groups[len(age_thresholds)]  # , age_group_cuis[len(age_thresholds)]


def find_sex(demographic):
    """Get sex/gender words out of the demographic."""
    sex_terms = {'male': ['male', 'man', 'boy'], 'female': ['female', 'woman', 'girl']}
    sex = demographic.split(' ')[3]
    return sex_terms[sex]
