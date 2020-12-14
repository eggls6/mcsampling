#############################################
#
# Query and transformation
# functions using JPL CNEOS API for asteroid
# orbital elements
#
# S. Eggl 20200123
############################################


import numpy as np
import json
import requests
import spiceypy as sp


def query_cometary_elements(tname):
    """Query JPL CNEOS API for cometary orbital elements of a minor planet at a
    given epoch from NASA JPL's web API

    Dependencies: json, requests

    Input: string, object name or designation (e.g. 'Eros')
    Output: [epoch (JD), cometary elements
                        [e,q(au),tp(JD),node(deg),peri(deg),inc(deg)]]

    For details see:
    https://ssd-api.jpl.nasa.gov/doc/sbdb.html

    """

    url = "https://ssd-api.jpl.nasa.gov/sbdb.api?sstr="+tname+"&cov=mat&phys-par=true&full-prec=true"

    r = requests.request("GET", url)
    ast = json.loads(r.text)

    # absolute magnitude H
    # h_v = ast['phys_par'][0]['value']

    # epoch of orbital elements at Tref [JD]
    epoch_jd = float(ast['orbit']['epoch'])

    # orbital elements
    elem = ast['orbit']['elements']

    hdr = []
    val = []
    for i in range(len(elem)):
        hdr.append(elem[i]['name'])
        val.append(elem[i]['value'])

    # cometary orbital elements: e,q[au],tp[MJD],node[deg],peri[deg],inc[deg]
    idx = [0, 2, 7, 4, 5, 3]
    elements = []
    for i in idx:
        if(elem[i]['name'] == 'tp'):
            elements.append(float(elem[i]['value']))
        else:
            elements.append(float(elem[i]['value']))

    return [epoch_jd, elements]


def query_cometary_ele_and_cov(tname):
    """Query JPL CNEOS API for cometary orbital elements and
    the corresponding covariance matrix of a minor planet at a given epoch
    from NASA JPL's web API

    Dependencies: json, requests
    Input: string, object name or designation (e.g. 'Eros')
    Output: [epoch (JD), cometary elements [e,q(au),tp(JD),node(deg),peri(deg),
             inc(deg)],[6x6 or 9x9 covariance matrix]]

    For details see:
    https://ssd-api.jpl.nasa.gov/doc/sbdb.html
    """

    url = "https://ssd-api.jpl.nasa.gov/sbdb.api?sstr="+tname+"&cov=mat&phys-par=true&full-prec=true"

    r = requests.request("GET",url)
    ast = json.loads(r.text)

    # absolute magnitude H
    # h_v = ast['phys_par'][0]['value']

    # epoch of orbital elements at Tref [JD]
    epoch_jd = float(ast['orbit']['covariance']['epoch'])

    elem = ast['orbit']['covariance']['elements']
    hdr = []
    val = []
    for i in range(len(elem)):
        hdr.append(elem[i]['name'])
        val.append(float(elem[i]['value']))
    elements = val

    # square root covariance matrix for cometary orbital elements
    mat = (np.array(ast['orbit']['covariance']['data'])).astype(float)
    # print(mat)
    return [epoch_jd, elements, mat]


def query_kepler_elements(tname):
    """Query JPL CNEOS API for Keplerian orbital elements of a minor planet
    at a given epoch from NASA JPL's web API

    Dependencies: json, requests
    Input: string, object name or designation (e.g. 'Eros')
    Output: [epoch (JD), cometary elements
                         [a(au),e,i(deg),peri(deg),node(deg),M(deg)]]

    For details see:
    https://ssd-api.jpl.nasa.gov/doc/sbdb.html

    """

    url = "https://ssd-api.jpl.nasa.gov/sbdb.api?sstr="+tname+"&cov=mat&phys-par=true&full-prec=true"

    r = requests.request("GET", url)
    ast = json.loads(r.text)

    # absolute magnitude H
    # h_v = ast['phys_par'][0]['value']

    # epoch of orbital elements at Tref [JD]
    epoch_jd = float(ast['orbit']['epoch'])

    # orbital elements ['e', 'a', 'q', 'i', 'om', 'w', 'ma', 'tp', 'per', 'n', 'ad']
    elem = ast['orbit']['elements']

    hdr = []
    val = []
    for i in range(len(elem)):
        hdr.append(elem[i]['name'])
        val.append(elem[i]['value'])

    # cometary orbital elements: e,q[au],tp[MJD],node[deg],peri[deg],inc[deg]
    idx = [1, 0, 3, 5, 4, 6]
    print(hdr)
    elements = []
    for i in idx:
        if(elem[i]['name'] == 'tp'):
            elements.append(float(elem[i]['value']))
        else:
            elements.append(float(elem[i]['value']))

    return [epoch_jd, elements]


def cometary2keplerian(epoch, ele, mu=0.01720209895**2):
    """Convert cometary orbital elements to Keplerian orbital elements

    Parameters:
    -----------
    epoch ... epoch of cometary elements [JD]
    ele   ... cometary elements [e,q[au],tp[JD],node[deg],peri[deg],inc[deg]]

    Optional:
    ---------
    mu ... Gravitational parameter (e.g. k**2*(M+m))

    Returns:
    --------
    kep ... Keplerian orbital elements
            [a, e, i[deg], w[deg], node[deg], M[deg]]

    """
    pix2 = 2*np.pi
    a = ele[1]/(1.-ele[0])
    M = np.sqrt(mu/a**3)*(epoch-ele[2])
    while(M < 0):
        M = M+pix2
    while(M > pix2):
        M = M-pix2

    # a, e, i, w, node, M
    kep = [a, ele[0], ele[5], ele[4], ele[3], np.rad2deg(M)]

    return kep


def cometary2cartesian(epoch, com, mu=0.01720209895**2):
    """Uses spiceypy to convert cometary orbital elements to
    HELIOCENTRIC (!) Cartesian states

    Parameters:
    -----------
    epoch ... epoch of orbital elements [JD]
    com ... cometary element array [e,q,tp(JD),node(deg),peri(deg),inc(deg)]]
    mu ... Gravitational parameter

    Returns:
    --------
    cart ... Cartesian state (x,y,z,vx,vy,vz).

    External dependencies:
    ----------------------
    spiceypy (imported as sp)
    numpy (imported as np)
    cometary2keplerian

    """
    kep = cometary2keplerian(epoch, com, mu)

#   Input for spiceypy.conics:
#   q    = pericenter distance
#   e    = eccentricity
#   i    = inclination (deg)
#   node = longitude of the ascending node (deg)
#   w    = argument of pericenter (deg)
#   M    = mean anomaly at epoch (deg)
#   T0   = epoch
#   mu   = gravitational parameter

    cart = sp.conics(np.array([com[1], com[0], np.deg2rad(com[5]),
                     np.deg2rad(com[3]), np.deg2rad(com[4]),
                     np.deg2rad(kep[5]), 0, mu]), 0)
    return cart


def cartesian2keplerian(epoch, state, mu=0.01720209895**2):
    """Uses spiceypy to convert Cartesian states to Keplerian orbital elements

    Parameters:
    -----------
    epoch ... epoch of orbital elements [JD]
    state ... Cartesian state (x,y,z,vx,vy,vz)
    mu ... Gravitational parameter

    Returns:
    --------
    kep ... orbital elements array
            a    = pericenter distance
            e    = eccentricity
            i    = inclination (deg)
            w    = argument of pericenter (deg)
            node = longitude of the ascending node (deg)
            M    = mean anomaly at epoch (deg)
            T0   = epoch
            mu   = gravitational parameter

    External dependencies:
    ----------------------
    spiceypy (imported as sp)
    numpy (imported as np)

    """

#   Output for spiceypy.oscelt:
#   q    = pericenter distance
#   e    = eccentricity
#   i    = inclination (deg)
#   node = longitude of the ascending node (deg)
#   w    = argument of pericenter (deg)
#   M    = mean anomaly at epoch (deg)
#   T0   = epoch
#   mu   = gravitational parameter

    oscelt = sp.oscelt(state, epoch, mu)

    kep = []
    # semimajor axis a from q
    kep.append(oscelt[0]/(1-oscelt[1]))
    # eccentricity
    kep.append(oscelt[1])
    # inclination
    kep.append(np.rad2deg(oscelt[2]))
    # w: argument of pericenter
    kep.append(np.rad2deg(oscelt[4]))
    # node
    kep.append(np.rad2deg(oscelt[3]))
    # mean anomaly
    kep.append(np.rad2deg(oscelt[5]))

    return kep, epoch, mu


def keplerian2cartesian(epoch, kep, mu=0.01720209895**2):
    """Uses spiceypy to convert Keplerian orbital elements
    to HELIOCENTRIC (!) Cartesian states

    Parameters:
    -----------
    epoch ... epoch of orbital elements [JD]
    kep ... orbital element array [a,e,inc(deg),peri/w(deg),node(deg),M(deg)]
    mu ... Gravitational parameter

    Returns:
    --------
    cart ... Cartesian state (x,y,z,vx,vy,vz).

    External dependencies:
    ----------------------
    spiceypy (imported as sp)
    numpy (imported as np)


    """
    q = kep[0]*(1-kep[1])

#   Input for spiceypy.conics:
#   q    = pericenter distance
#   e    = eccentricity
#   i    = inclination (deg)
#   node = longitude of the ascending node (deg)
#   w    = argument of pericenter (deg)
#   M    = mean anomaly at epoch (deg)
#   T0   = epoch
#   mu   = gravitational parameter

    cart = sp.conics(np.array([q, kep[1], np.deg2rad(kep[2]),
                               np.deg2rad(kep[4]), np.deg2rad(kep[3]),
                               np.deg2rad(kep[5]), 0, mu]), 0)
    return cart
