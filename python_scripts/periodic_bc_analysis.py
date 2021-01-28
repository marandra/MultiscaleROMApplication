# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
    StructuralMechanicsAnalysis,
)


class PBCAnalysis(StructuralMechanicsAnalysis):
    def __init__(self, model, project_parameters):

        self.boundary_mp_name = project_parameters["rve_settings"][
            "boundary_mp_name"
        ].GetString()
        self.averaging_mp_name = project_parameters["rve_settings"][
            "averaging_mp_name"
        ].GetString()

        super(PBCAnalysis, self).__init__(model, project_parameters)

    def ModifyInitialGeometry(self):
        # Here populate the submodelparts to be used for periodicity
        super(PBCAnalysis, self).ModifyInitialGeometry()

        boundary_mp = self.model[self.boundary_mp_name]
        averaging_mp = self.model[self.averaging_mp_name]

        # Construct auxiliary modelparts
        self.min_corner, self.max_corner = self._DetectBoundingBox(averaging_mp)
        self._ConstructFaceModelParts(self.min_corner, self.max_corner, boundary_mp)

    def Initialize(self):
        # construct MPCs according to the provided strain
        super(PBCAnalysis, self).Initialize()
        boundary_mp = self.model[self.boundary_mp_name]
        averaging_mp = self.model[self.averaging_mp_name]
        strain = KratosMultiphysics.Matrix(3, 3, 0.0)
        self._ApplyPeriodicity(strain, averaging_mp, boundary_mp)

    def _DetectBoundingBox(self, mp):
        min_corner = KratosMultiphysics.Array3()
        min_corner[0] = 1e20
        min_corner[1] = 1e20
        min_corner[2] = 1e20

        max_corner = KratosMultiphysics.Array3()
        max_corner[0] = -1e20
        max_corner[1] = -1e20
        max_corner[2] = -1e20

        for node in mp.Nodes:
            x = node.X
            min_corner[0] = min(min_corner[0], x)
            max_corner[0] = max(max_corner[0], x)

            y = node.Y
            min_corner[1] = min(min_corner[1], y)
            max_corner[1] = max(max_corner[1], y)

            z = node.Z
            min_corner[2] = min(min_corner[2], z)
            max_corner[2] = max(max_corner[2], z)

        KratosMultiphysics.Logger.PrintInfo("Periodic BC", "Boundng box detected")
        KratosMultiphysics.Logger.PrintInfo("Periodic BC", "Min. corner = ", min_corner)
        KratosMultiphysics.Logger.PrintInfo("Periodic BC", "Max. corner = ", max_corner)

        return min_corner, max_corner

    def __PopulateMp(self, face_name, coordinate, component, eps, mp):
        if mp.NumberOfConditions() == 0:
            raise Exception("Boundary_mp is expected to have conditions and has none")

        mp = mp.GetRootModelPart()

        if not mp.HasSubModelPart(face_name):
            mp.CreateSubModelPart(face_name)
        face_mp = mp.GetSubModelPart(face_name)

        for cond in mp.Conditions:
            xc = cond.GetGeometry().Center()
            if abs(xc[component] - coordinate) < eps:
                face_mp.AddCondition(cond)

        node_ids = set()
        for cond in face_mp.Conditions:
            for node in cond.GetNodes():
                if not node.Is(KratosMultiphysics.SLAVE):
                    node_ids.add(node.Id)
                    node.Set(KratosMultiphysics.SLAVE)

        face_mp.AddNodes(list(node_ids))
        return face_mp

    def _ConstructFaceModelParts(self, min_corner, max_corner, mp):

        eps = 0.0001 * (max_corner[0] - min_corner[0]) / mp.NumberOfNodes()

        KratosMultiphysics.VariableUtils().SetFlag(
            KratosMultiphysics.SLAVE, False, mp.Nodes
        )
        KratosMultiphysics.VariableUtils().SetFlag(
            KratosMultiphysics.MASTER, False, mp.Nodes
        )

        # Populate the slave faces
        self.max_x_face = self.__PopulateMp("max_x_face", max_corner[0], 0, eps, mp)
        self.max_y_face = self.__PopulateMp("max_y_face", max_corner[1], 1, eps, mp)
        self.max_z_face = self.__PopulateMp("max_z_face", max_corner[2], 2, eps, mp)

        # First populate the master faces (min)
        self.min_x_face = self.__PopulateMp("min_x_face", min_corner[0], 0, eps, mp)
        self.min_y_face = self.__PopulateMp("min_y_face", min_corner[1], 1, eps, mp)
        self.min_z_face = self.__PopulateMp("min_z_face", min_corner[2], 2, eps, mp)

        if self.min_x_face.NumberOfConditions() == 0:
            raise Exception("min_x_face has 0 conditions")
        if self.min_y_face.NumberOfConditions() == 0:
            raise Exception("min_y_face has 0 conditions")
        if self.min_z_face.NumberOfConditions() == 0:
            raise Exception("min_z_face has 0 conditions")

    def _ApplyPeriodicity(self, strain, volume_mp, boundary_mp):
        # clear
        for constraint in volume_mp.GetRootModelPart().MasterSlaveConstraints:
            constraint.Set(KratosMultiphysics.TO_ERASE)
        volume_mp.GetRootModelPart().RemoveMasterSlaveConstraintsFromAllLevels(
            KratosMultiphysics.TO_ERASE
        )

        dx = self.max_corner[0] - self.min_corner[0]
        dy = self.max_corner[1] - self.min_corner[1]
        dz = self.max_corner[2] - self.min_corner[2]

        periodicity_utility = KratosMultiphysics.StructuralMechanicsApplication.RVEPeriodicityUtility(
            self._GetSolver().GetComputingModelPart()
        )

        # assign periodicity to faces
        periodicity_utility.AssignPeriodicity(
            self.min_x_face,
            self.max_x_face,
            strain,
            KratosMultiphysics.Vector([dx, 0.0, 0.0]),
            0.001
        )
        #periodicity_utility.AssignPeriodicity(
        #    self.min_y_face,
        #    self.max_y_face,
        #    strain,
        #    KratosMultiphysics.Vector([0.0, dy, 0.0]),
        #    0.001
        #)
        periodicity_utility.AssignPeriodicity(
            self.min_z_face,
            self.max_z_face,
            strain,
            KratosMultiphysics.Vector([0.0, 0.0, dz]),
            0.001
        )

        periodicity_utility.Finalize(KratosMultiphysics.DISPLACEMENT)

        # start from the exact solution in the case of a constant strain
        x = KratosMultiphysics.Array3()
        for node in volume_mp.Nodes:
            x[0] = node.X0
            x[1] = node.Y0
            x[2] = node.Z0
            d = strain * x
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, d)
