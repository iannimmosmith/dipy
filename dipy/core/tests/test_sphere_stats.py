import numpy as np
import numpy.testing as nt
from dipy.data import get_sphere


from dipy.core.sphere_stats import (evaluate_sphere_eqzone, 
                                    evaluate_sphere_coarseness)

def test_evaluate_eqzone():
    vertices, faces = get_sphere('symmetric362')
    eqc, polc = evaluate_sphere_eqzone(vertices, width=5)
    nt.assert_almost_equal(np.mean(eqc), 30.165745856353592)
    
    #vertices2, faces2 = get_sphere('symmetric724')
    #eqc2, polc2 = evaluate_sphere_eqzone(vertices2, width=5)
    
    coarse = evaluate_sphere_coarseness(vertices, faces)
    #coarse2 = evaluate_sphere_coarseness(vertices2, faces2)
    nt.assert_almost_equal(coarse, 0.43842554416189522)
    
if __name__ == "__main__":
    nt.run_module_suite()
