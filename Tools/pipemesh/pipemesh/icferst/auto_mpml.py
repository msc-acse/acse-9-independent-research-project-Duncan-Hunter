import xml.etree.ElementTree as ET
import copy
import os
loc = os.path.dirname(os.path.abspath(__file__))

class AutoMPML():
    def __init__(self, fname=loc+"/3d_pipe_FEM.mpml"):
        self.mpml_tree = ET.parse(fname)
        self.mpml_root = self.mpml_tree.getroot()

        self.sim_name = self.mpml_root[0]
        self.msh_options = self.mpml_root[2]
        self.io_options = self.mpml_root[4]
        self.timestepping = self.mpml_root[5]
        self.material_phase = self.mpml_root[6]
        self.mesh_adaptivity = self.mpml_root[7]

    def set_sim_name(self, filename):
        self.sim_name[0].text = filename

    def set_msh_options(self, msh_file):
        mesh = self.msh_options[1]
        mesh[0].attrib['file_name'] = msh_file

    def set_io_options(self, dump_period=0.1):
        if dump_period < 0.1:
            print("Warning: dump period is less than maximum timestep.")
        dump_period = str(dump_period)
        dump_p = self.io_options[1][0][0]
        dump_p.text = str(dump_period)

    def set_timestepping(self, finish_time, timestep=0.005, CFL_no=2):
        tstep = self.timestepping[1][0]
        tstep.text = str(timestep)
        ftime = self.timestepping[2][0]
        ftime.text = str(finish_time)

        if CFL_no > 4:
            raise ValueError("CFL number too high")
        elif CFL_no > 2:
            print("Warning: High CFL number")
        elif CFL_no <= 0:
            raise ValueError("CFL number too low")

        cfl = self.timestepping[3][0][0]
        cfl.text = str(CFL_no)

    def set_material_properties(self, density=1e3, viscosity=1e-3):
        density_o = self.material_phase[0][0][0][0]
        density_o.text = str(density)
        viscosity_o = self.material_phase[0][1][0][0][0][0][0][0]
        viscosity_o.text = str(viscosity)
        
    def set_inlets(self, phys_ids, velocities, directions=None):
        """Sets the properties for inlets.
        
        Args:
            phys_ids: (list of ints) physical ids of inlet surfaces.
            velocities: (list of xyz vector lists) velocities
                for respective physical ids. E.g. velocities[0] corresponds to 
                phys_ids[0]. If only one vector given, 
            directions: (list of xyz vector lists) outward
                directions of surfaces. Used to check if direction aligns with velocity.
                Reversed for inlets. Can be left.
        """
        def set_velo(i1_elem, i1_mom_elem, i1_ad_elem, phys_id, velocity):
            """eleme is boundary_conditions, name=inlet."""
            def set_vec_comp(elem, velo):
                for i in range(3):
                    elem[i][0][0].text = str(velo[i])
            i1_elem.attrib['name'] = "inlet_{}".format(phys_id)
            i1_elem[0][0].text = str(phys_id)
            set_vec_comp(i1_elem[1][1], velocity)      
            
            i1_mom_elem.attrib['name'] = "inlet_{}_mom".format(phys_id)
            i1_mom_elem[0][0].text = str(phys_id)
            set_vec_comp(i1_mom_elem[1][0], velocity)

            i1_ad_elem.attrib['name'] = "inlet_{}_ad".format(phys_id)
            i1_ad_elem[0][0].text = str(phys_id)
            set_vec_comp(i1_ad_elem[1][1], velocity)

        velo = self.material_phase[2]

        n_inlets = len(phys_ids)
        set_velo(velo[0][2], velo[0][3], velo[0][4], phys_ids[0], velocities[0])
        if n_inlets > 1:
            if len(velocities) < n_inlets:
                print("Not enough velocities given. Using first velocity for all inlets.")
                velocities = [velocities[0]] * n_inlets
            for phys_id, velocity in zip(phys_ids[1:], velocities[1:]):
                i1_copy = copy.deepcopy(velo[0][2])
                i1_mom_copy = copy.deepcopy(velo[0][3])
                i1_ad_copy = copy.deepcopy(velo[0][4])
                set_velo(i1_copy, i1_mom_copy, i1_ad_copy, phys_id, velocity)
                velo[0].append(i1_copy)
                velo[0].append(i1_mom_copy)
                velo[0].append(i1_ad_copy)
    
    def set_outlets(self, phys_ids):
        pressure = self.material_phase[1]
        out_ids = pressure[0][2][0][0]
        out_ids.attrib['shape'] = len(phys_ids)
        text = str(phys_ids[0])
        for i in phys_ids[1:]:
            text += " {}".format(int(i))
        out_ids.text = text
    
    def set_no_slip(self, phys_ids):
        cyl_mom = self.material_phase[2][0][5]
        cyl_ns_visc = self.material_phase[2][0][6]
        cyl = self.material_phase[2][0][7]
        text = str(phys_ids[0])
        for i in phys_ids[1:]:
            text += " {}".format(int(i))

        for cond in [cyl_mom, cyl_ns_visc, cyl]:
            out_ids = cond[0][0]
            out_ids.attrib['shape'] = len(phys_ids)
            out_ids.text = text
    
    def set_mesh_adaptivity(self,
                            min_size=0.01,
                            max_size=0.5,
                            max_no_nodes=300000,
                            t_adapt_delay=0.5,
                            aspect_ratio=5):
        mnn = self.mesh_adaptivity[0][1][0]
        mnn.text = str(int(max_no_nodes))
        
        def ani_sym_matrix_text(value):
            text = "{} 0.0 0.0 0.0 {} 0.0 0.0 0.0 {}".format(value, value, value)
            return text
        min_s = self.mesh_adaptivity[0][3][0][0][0]
        min_s.text = ani_sym_matrix_text(min_size)
        max_s = self.mesh_adaptivity[0][4][0][0][0]
        max_s.text = ani_sym_matrix_text(max_size)
        asp_ratio = self.mesh_adaptivity[0][4][0]
        asp_ratio.text = str(aspect_ratio)
        t_a_d = self.mesh_adaptivity[0][5][0]
        t_a_d.text = str(t_adapt_delay)
    
    def write_mpml(self):
        self.mpml_tree.write("test.mpml", 'utf-8', xml_declaration=True)
