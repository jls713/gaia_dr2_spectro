from spectro_data.apogee import load_and_match as load_apogee
from spectro_data.test_apogeetgas import load_and_match as load_apogeetgas
from spectro_data.galah import load_and_match as load_galah
from spectro_data.ges import load_and_match as load_ges
from spectro_data.lamost import load_and_match as load_lamost
from spectro_data.rave import load_and_match as load_rave
from distances import *

DR1 = True
TEST = False # True
TESTNUMBER = 100
npool = 10
random_seed = 52

with open('config.json') as config:
    output_folder = json.load(config)['dir']['output_folder']

# RAVE


def run_rave_on():
    rave = load_rave(output_folder + 'RAVEON_input.hdf5', use_dr1=DR1)
    if TEST:
        rave = rave.sample(n=TESTNUMBER,
                           random_state=random_seed).reset_index(drop=True)
    run_distance_pipeline(rave,
                          output_folder + 'RAVE_DR5_distances.hdf5',
                          'raveid', 'RAVE DR5',
                          npool=npool)


def run_rave_dr5():
    rave = load_rave(use_dr1=DR1)
    if TEST:
        rave = rave.sample(n=TESTNUMBER,
                           random_state=random_seed).reset_index(drop=True)
    run_distance_pipeline(rave,
                          output_folder + 'RAVE_Cannon_distances.hdf5',
                          'raveid', 'RAVEon',
                          npool=npool)

# Gaia-ESO


def run_ges():
    ges = load_ges(use_dr1=DR1)
    if TEST:
        ges = ges.sample(n=TESTNUMBER,
                         random_state=random_seed).reset_index(drop=True)
    run_distance_pipeline(ges,
                          output_folder + 'GES_distances.hdf5',
                          'CNAME', 'GES',
                          npool=npool)

# APOGEE


def run_apogee():
    apogee = load_apogee(use_dr1=DR1)
    if TEST:
        apogee = apogee.sample(n=TESTNUMBER,
                               random_state=random_seed).reset_index(drop=True)
    run_distance_pipeline(apogee,
                          output_folder + 'APOGEE_distances.hdf5',
                          'APOGEE_ID', 'APOGEE',
                          npool=npool)

# LAMOST


def run_lamost():
    lamost = load_lamost(use_dr1=DR1)
    if TEST:
        lamost = lamost.sample(n=TESTNUMBER,
                               random_state=random_seed).reset_index(drop=True)
    run_distance_pipeline(lamost,
                          output_folder + 'LAMOST_distances.hdf5',
                          'obsid', 'LAMOST',
                          npool=npool)

# GALAH


def run_galah():
    galah = load_galah(use_dr1=DR1)
    if TEST:
        galah = galah.sample(n=TESTNUMBER,
                             random_state=random_seed).reset_index(drop=True)
    run_distance_pipeline(galah,
                          output_folder + 'GALAH_distances.hdf5',
                          'sobject_id', 'GALAH',
                          npool=npool)

# APOGEE TGAS TEST


def run_apogeetgas_test():
    apogee = load_apogeetgas()
    print len(apogee)
    if TEST:
        apogee = apogee.sample(n=TESTNUMBER,
                               random_state=random_seed).reset_index(drop=True)
    print len(apogee)
    run_distance_pipeline(apogee,
                          output_folder +
                          'APOGEETGAS_TEST_distances_random100.hdf5',
                          'location_id', 'APOGEE TGAS',
                          npool=npool)


def run_apogeetgas_test_withGaia():
    apogee = load_apogeetgas()
    apogee['eG'] = 0.025
    apogee['mag_use'] = \
        apogee.applymap(
            lambda x: np.array(['J', 'H', 'K', 'G']))['mag_use']
    if TEST:
        apogee = apogee.sample(n=TESTNUMBER,
                               random_state=random_seed).reset_index(drop=True)
    run_distance_pipeline(apogee,
                          output_folder +
                          'APOGEETGAS_TEST_distances_random100_G.hdf5',
                          'location_id', 'APOGEE TGAS',
                          npool=npool)


if __name__ == '__main__':
    # run_rave_on()
    # run_rave_dr5()
    # run_ges()
    # run_apogee()
    # run_lamost()
    run_galah()
    # run_apogeetgas_test()
    # run_apogeetgas_test_withGaia()
