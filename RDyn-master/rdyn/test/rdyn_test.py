import unittest
import shutil
from rdyn.alg.RDyn_v2 import RDynV2


class RDynTestCase(unittest.TestCase):

    def test_rdyn_simplified(self):

        print("1")
        rdb = RDynV2(size=500, iterations=100)
        rdb.execute(simplified=True)
        print("2")
        rdb = RDynV2(size=500, iterations=100, max_evts=2)
        rdb.execute(simplified=True)
        print("3")
        rdb = RDynV2(size=500, iterations=100, new_node=0.1, del_node=0.1, max_evts=2, paction=0.8)
        rdb.execute(simplified=False)
        print("Done")

        shutil.rmtree("results")


if __name__ == '__main__':
    unittest.main()
