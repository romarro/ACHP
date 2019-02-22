import numpy as np
from ACHP.PHEHX import CP, PHEHXClass

def test_pf_BuildEnthalpyLists_with_water_and_R1234yf_should_return_two_zone():
    """Test to see if the entalhyList for parallel flow returns the expected values for a specific case"""

    # arrange     
    ref = CP.AbstractState('HEOS','R1234yf')
    
    ref.update(CP.PQ_INPUTS, 1.2e5, 0.4)
    ref_inlet = {
        'AS_c':ref,
        'mdot_c':0.05,
        'pin_c':ref.p(),
        'hin_c':ref.hmass()
    }

    wg = CP.AbstractState('INCOMP','Water')
    wg.update(CP.PT_INPUTS, 2e5, 40.0 + 273.15)
    wg_inlet = {
        'AS_h':wg,
        'mdot_h':0.25 * wg.rhomass() * 1e-3, # din l/s -> kg/s
        'pin_h':2e5,
        'hin_h':wg.hmass()
    }
    
    params = {}
    params.update(**ref_inlet, **wg_inlet)
    phx = PHEHXClass(**params)

    Qdot = 5553.217 # W
    ref_hlist_expected = [239402.32449262452, 346315.0179186873, 350466.66287818813]
    wg_hlist_expected = [83531.24578146616, 62019.66312319436, 61184.32308215549]
    # act 

    ref_hlist_actual, wg_hlist_actual = phx._pf_BuildEnthalpyLists(Qdot)

    # assert
    assert len(ref_hlist_expected) == len(ref_hlist_actual), 'refrigerant list length not ok'
    assert len(wg_hlist_expected) == len(wg_hlist_actual), 'water list length not ok'
    assert np.allclose(ref_hlist_expected, ref_hlist_actual,atol=1e-3), 'refrigerant list not ok'
    assert np.allclose(wg_hlist_expected, wg_hlist_actual, atol=1e-3), 'water list not ok'


def test_PHEHX_Parallel_flow_with_water_and_r1234yf():

    # arrange
    ref = CP.AbstractState('HEOS','R1234yf')
    
    ref.update(CP.PQ_INPUTS, 3.8e5, 0.4)
    ref_inlet = {
        'AS_c':ref,
        'mdot_c':0.05,
        'pin_c':ref.p(),
        'hin_c':ref.hmass()
    }

    wg = CP.AbstractState('INCOMP','Water')
    wg.update(CP.PT_INPUTS, 2e5, 40.0 + 273.15)
    wg_inlet = {
        'AS_h':wg,
        'mdot_h':0.25 * wg.rhomass() * 1e-3, # din l/s -> kg/s
        'pin_h':2e5,
        'hin_h':wg.hmass()
    }
    
    #Geometric parameters
    params = {
        #Geometric parameters
        'Bp' : 0.101,
        'Lp' : 0.455, #Center-to-center distance between ports
        'Nplates' : 46,
        'PlateAmplitude' : 0.00102, #[m]
        'PlateThickness' : 0.0003, #[m]
        'PlateWavelength' : 0.00626, #[m]
        'InclinationAngle' : 65/180*np.pi,#[rad]
        'PlateConductivity' : 15.0, #[W/m-K]
        'Rp': 1.0, #[microns] Surface roughness
        'MoreChannels' : 'Hot', #Which stream gets the extra channel, 'Hot' or 'Cold'
        'FlowConfig' : 'Parallel', # added to enable parallel flow / or Counter flow
    
        'Verbosity':0,
        
        'h_tp_cold_tuning':1,
        'h_tp_hot_tuning':0.5,
        'DP_hot_tuning':1,
        'DP_cold_tuning':1
    }
    params.update(**ref_inlet, **wg_inlet)
    phx = PHEHXClass(**params)
    # act

    phx.Calculate()

    
    # assert
    assert True