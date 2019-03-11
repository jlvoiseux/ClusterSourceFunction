var cluster;
var cluster_length;
var curr_atom;
var curr_x;
var curr_y;
var curr_z;
var curr_atom_bond;
var curr_x_bond;
var curr_y_bond;
var curr_z_bond;
var Dx;
var Dy;
var Dz;
var rotx;
var roty;
var rotz;
var max_dist_bond = 3;
var bond_radius = 0.1;
var scale_const = 1;
var size_const = 0.3;
var size_eq = 5;
var min_x = Infinity;
var min_y = Infinity;
var min_z = Infinity;
var bond_array;
var controls;
var mouse = new THREE.Vector2(),
    INTERSECTED;
var CPraw;
var CPproperties;
var click_array;
var atom_id;
var gui;
var cp_nam_item_fun = function () {
    this.Atom_Name = "0";
};
var cp_pos_item_fun = function () {
    this.Atom_Position_Bohr = "0";
};
var cp_den_item_fun = function () {
    this.Electron_Density = "0";
};
var cp_lag_item_fun = function () {
    this.Lagrangian_Kinetic_Energy = "0";
};
var cp_ham_item_fun = function () {
    this.Hamiltonian_Kinetic_Energy = "0";
};
var cp_pot_item_fun = function () {
    this.Potential_Energy_Density = "0";
};
var cp_end_item_fun = function () {
    this.Energy_Density = "0";
};
var cp_lap_item_fun = function () {
    this.Laplacian_Of_Electron_Density = "0";
};
var cp_loc_item_fun = function () {
    this.Electron_Localization_Function = "0";
};
var source_item_fun = function () {
    this.Source_Contribution = "0";
};
var bond_item_fun = function () {
    this.Bond_Distance = max_dist_bond;
};
var source_index_item_fun = function () {
    this.IP_Index = 1;
};
var source_name_item_fun = function () {
    this.IP_Name = "0";
};

var cp_nam_item;
var cp_pos_item;
var cp_den_item;
var cp_lag_item;
var cp_ham_item;
var cp_pot_item;
var cp_end_item;
var cp_lap_item;
var cp_loc_item;
var source_item;
var bond_item;
var source_index_item;
var source_name_item;

var source_raw;
var source_contrib;
var source_points;

var prev_index;



var clock = new THREE.Clock();

var atom_geometry = new THREE.SphereGeometry(1, 125, 125);
var atom_material = new THREE.MeshBasicMaterial({
    color: 0x00ff00
});
var bond_geometry = new THREE.CylinderGeometry(bond_radius, bond_radius, 1, 125);
var bond_material = new THREE.MeshBasicMaterial({
    color: 0x00ff00
});
var scene = new THREE.Scene();
scene.background = new THREE.Color(0xf0f0f0);
var camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
camera.up = new THREE.Vector3(0, 1, 0);
var renderer = new THREE.WebGLRenderer();

var light = new THREE.DirectionalLight(0xffffff, 1);
light.position.set(1, 1, 1).normalize();

var light2 = new THREE.DirectionalLight(0xffffff, 1);
light2.position.set(-1, -1, -1).normalize();

raycaster = new THREE.Raycaster();

renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);
document.addEventListener('mousemove', onDocumentMouseMove, false);
controls = new THREE.TrackballControls(camera, renderer.domElement);


load_cluster();

camera.position.z = 5;

function load_cluster() {
    var loader = new THREE.FileLoader();

    //load a text file and output the result to the console
    loader.load(
        // resource URL
        'cluster_viz.txt',

        // onLoad callback
        function (data) {
            // output the text to the console

            cluster = data.split("\n")
            cluster_length = cluster.length;            
            animate();
        },

        // onProgress callback
        function (xhr) {
            console.log((xhr.loaded / xhr.total * 100) + '% loaded');
        },

        // onError callback
        function (err) {
            console.error('An error happened');
        }
    );

    loader.load(
        // resource URL
        'CPprop.txt',

        // onLoad callback
        function (data) {
            // output the text to the console
            CPraw = data.split("\n")
        },

        // onProgress callback
        function (xhr) {
            console.log((xhr.loaded / xhr.total * 100) + '% loaded');
        },

        // onError callback
        function (err) {
            console.error('An error happened');
        }
    );

    loader.load(
        // resource URL
        'source.csv',

        // onLoad callback
        function (data) {
            // output the text to the console
            source_raw = data.split("\n")           
            make_molecule_init();
            make_CP_prop();
            make_GUI();
            

        },

        // onProgress callback
        function (xhr) {
            console.log((xhr.loaded / xhr.total * 100) + '% loaded');
        },

        // onError callback
        function (err) {
            console.error('An error happened');
        }
    );

}



function animate() {
    requestAnimationFrame(animate);
    var delta = clock.getDelta();
    if(max_dist_bond !=3 && bond_item.Bond_Distance != max_dist_bond){
        while(scene.children.length > 0){ 
            scene.remove(scene.children[0]); 
        }
        make_molecule();
    }    
    controls.update(delta);
    render();
}

function render() {
    raycaster.setFromCamera(mouse, camera);
    var intersects = raycaster.intersectObjects(scene.children);
    click_array = new Array();
    if (intersects.length > 0) {
        if (INTERSECTED != intersects[0].object) {
            if (INTERSECTED) INTERSECTED.material.emissive.setHex(INTERSECTED.currentHex);
            INTERSECTED = intersects[0].object;
            for (var i = 0; i < CPproperties.length; i++) {
                if (parseInt(INTERSECTED.userData.id) == (CPproperties[i][0])) {
                    update_GUI(CPproperties[i][1], CPproperties[i][2], CPproperties[i][3], CPproperties[i][4], CPproperties[i][5], CPproperties[i][6], CPproperties[i][7], CPproperties[i][8], CPproperties[i][9], CPproperties[i][9 + source_index_item.IP_Index], source_points[source_index_item.IP_Index-1]);
                }
            }
            INTERSECTED.currentHex = INTERSECTED.material.emissive.getHex();
            INTERSECTED.material.emissive.setHex(0xff0000);
        }
    } else {
        if (INTERSECTED) INTERSECTED.material.emissive.setHex(INTERSECTED.currentHex);
        INTERSECTED = null;
    }
    renderer.render(scene, camera);
}

function make_molecule_init() {

    scene.add(light);
    scene.add(light2);
    
    for (var i = 0; i < cluster_length; i++) {
        curr_atom = cluster[i].split(/\s+/);
        curr_x = parseFloat(curr_atom[2]) * scale_const;
        curr_y = parseFloat(curr_atom[3]) * scale_const;
        curr_z = parseFloat(curr_atom[4]) * scale_const;

        if (Math.abs(curr_x) < min_x) {
            min_x = curr_x;
        }
        if (Math.abs(curr_y) < min_y) {
            min_y = curr_y;
        }
        if (Math.abs(curr_z) < min_z) {
            min_z = curr_z;
        }
    }


    bond_array = new Array(cluster_length);
    for (var i = 0; i < cluster_length; i++) {
        bond_array[i] = new Array(cluster_length);
    }

    for (var i = 0; i < cluster_length; i++) {
        curr_atom = cluster[i].split(/\s+/);
        for (var j = 0; j < cluster_length; j++) {
            if (i != j && j > i) {
                curr_x = parseFloat(curr_atom[2]) * scale_const - 2 * min_x;
                curr_y = parseFloat(curr_atom[3]) * scale_const - 2 * min_y;
                curr_z = parseFloat(curr_atom[4]) * scale_const - 2 * min_z;
                curr_atom_bond = cluster[j].split(/\s+/);
                curr_x_bond = parseFloat(curr_atom_bond[2]) * scale_const - 2 * min_x;
                curr_y_bond = parseFloat(curr_atom_bond[3]) * scale_const - 2 * min_y;
                curr_z_bond = parseFloat(curr_atom_bond[4]) * scale_const - 2 * min_z;
                Dx = curr_x - curr_x_bond;
                Dy = curr_y - curr_y_bond;
                Dz = curr_z - curr_z_bond;
                rotx = Math.atan(Dz / Dy);
                if (Math.abs(Dy) < 0.005 && Math.abs(Dz) < 0.005) {
                    rotx = 0
                } else if (Math.abs(Dy) < 0.005 && Math.abs(Dx) < 0.005) {
                    rotx = Math.PI / 2;
                } else if (Math.abs(Dy) < 0.005 && Dy > 0) {
                    rotx = Math.PI / 2;
                } else if (Math.abs(Dy) < 0.005 && Dy < 0) {
                    rotx = 0;
                }
                roty = Math.atan(Dz / Dx);
                rotz = Math.atan(Dx / Dy);
                if (Math.abs(Dy) < 0.005 && Math.abs(Dz) < 0.005) {
                    rotz = Math.PI / 2;
                } else if (Math.abs(Dy) < 0.005 && Math.abs(Dx) < 0.005) {
                    rotz = 0;
                } else if (Math.abs(Dy) < 0.005) {
                    rotz = Math.PI / 4;
                }
                bond_array[i][j] = [Math.sqrt(Math.pow(Dx, 2) + Math.pow(Dy, 2) + Math.pow(Dz, 2)), rotz, rotx, roty, (curr_x + curr_x_bond) / 2, (curr_y + curr_y_bond) / 2, (curr_z + curr_z_bond) / 2, curr_x - curr_x_bond, curr_y - curr_y_bond, curr_z - curr_z_bond]
            } else {
                bond_array[i][j] = [0, 0, 0]
            }
        }
    }

    atom_id = new Array(cluster_length);

    for (var i = 0; i < cluster_length; i++) {
        curr_atom = cluster[i].split(/\s+/);

        curr_name = curr_atom[0];
        curr_load = parseFloat(curr_atom[1]);
        curr_x = parseFloat(curr_atom[2]) * scale_const - 2 * min_x;
        curr_y = parseFloat(curr_atom[3]) * scale_const - 2 * min_y;
        curr_z = parseFloat(curr_atom[4]) * scale_const - 2 * min_z;

        var atom = new THREE.Mesh(atom_geometry, new THREE.MeshLambertMaterial({
            color: curr_load / 100 * curr_load / 100 * 0xffffff
        }));

        atom.position.x = curr_x;
        atom.position.y = curr_y;
        atom.position.z = curr_z;
        atom.scale.x = Math.log(curr_load + size_eq) * size_const;
        atom.scale.y = Math.log(curr_load + size_eq) * size_const;
        atom.scale.z = Math.log(curr_load + size_eq) * size_const;
        atom.userData.id = (i + 1).toString();
        atom_id[i] = [atom.userData.id, parseFloat(curr_atom[2]), parseFloat(curr_atom[3]), parseFloat(curr_atom[4]), curr_atom[0]];
        scene.add(atom);
    }

    for (var i = 0; i < cluster_length; i++) {
        for (var j = 0; j < cluster_length; j++) {
            if (bond_array[i][j][0] != 0 && bond_array[i][j][0] < max_dist_bond) {
                var bond = new THREE.Mesh(bond_geometry, new THREE.MeshLambertMaterial({
                    color: 0.9 * 0xffffff
                }));
                bond.position.x = bond_array[i][j][4];
                bond.position.y = bond_array[i][j][5];
                bond.position.z = bond_array[i][j][6];
                bond.scale.y = bond_array[i][j][0];
                var vector = new THREE.Vector3(bond_array[i][j][7], bond_array[i][j][8], bond_array[i][j][9]);
                var axis = new THREE.Vector3(0, 1, 0);
                bond.quaternion.setFromUnitVectors(axis, vector.clone().normalize());
                scene.add(bond);

            }
        }
    }
}

function make_molecule() {

    scene.add(light);
    scene.add(light2);
    
    for (var i = 0; i < cluster_length; i++) {
        curr_atom = cluster[i].split(/\s+/);
        curr_x = parseFloat(curr_atom[2]) * scale_const;
        curr_y = parseFloat(curr_atom[3]) * scale_const;
        curr_z = parseFloat(curr_atom[4]) * scale_const;

        if (Math.abs(curr_x) < min_x) {
            min_x = curr_x;
        }
        if (Math.abs(curr_y) < min_y) {
            min_y = curr_y;
        }
        if (Math.abs(curr_z) < min_z) {
            min_z = curr_z;
        }
    }


    bond_array = new Array(cluster_length);
    for (var i = 0; i < cluster_length; i++) {
        bond_array[i] = new Array(cluster_length);
    }

    for (var i = 0; i < cluster_length; i++) {
        curr_atom = cluster[i].split(/\s+/);
        for (var j = 0; j < cluster_length; j++) {
            if (i != j && j > i) {
                curr_x = parseFloat(curr_atom[2]) * scale_const - 2 * min_x;
                curr_y = parseFloat(curr_atom[3]) * scale_const - 2 * min_y;
                curr_z = parseFloat(curr_atom[4]) * scale_const - 2 * min_z;
                curr_atom_bond = cluster[j].split(/\s+/);
                curr_x_bond = parseFloat(curr_atom_bond[2]) * scale_const - 2 * min_x;
                curr_y_bond = parseFloat(curr_atom_bond[3]) * scale_const - 2 * min_y;
                curr_z_bond = parseFloat(curr_atom_bond[4]) * scale_const - 2 * min_z;
                Dx = curr_x - curr_x_bond;
                Dy = curr_y - curr_y_bond;
                Dz = curr_z - curr_z_bond;
                rotx = Math.atan(Dz / Dy);
                if (Math.abs(Dy) < 0.005 && Math.abs(Dz) < 0.005) {
                    rotx = 0
                } else if (Math.abs(Dy) < 0.005 && Math.abs(Dx) < 0.005) {
                    rotx = Math.PI / 2;
                } else if (Math.abs(Dy) < 0.005 && Dy > 0) {
                    rotx = Math.PI / 2;
                } else if (Math.abs(Dy) < 0.005 && Dy < 0) {
                    rotx = 0;
                }
                roty = Math.atan(Dz / Dx);
                rotz = Math.atan(Dx / Dy);
                if (Math.abs(Dy) < 0.005 && Math.abs(Dz) < 0.005) {
                    rotz = Math.PI / 2;
                } else if (Math.abs(Dy) < 0.005 && Math.abs(Dx) < 0.005) {
                    rotz = 0;
                } else if (Math.abs(Dy) < 0.005) {
                    rotz = Math.PI / 4;
                }
                bond_array[i][j] = [Math.sqrt(Math.pow(Dx, 2) + Math.pow(Dy, 2) + Math.pow(Dz, 2)), rotz, rotx, roty, (curr_x + curr_x_bond) / 2, (curr_y + curr_y_bond) / 2, (curr_z + curr_z_bond) / 2, curr_x - curr_x_bond, curr_y - curr_y_bond, curr_z - curr_z_bond]
            } else {
                bond_array[i][j] = [0, 0, 0]
            }
        }
    }

    atom_id = new Array(cluster_length);

    for (var i = 0; i < cluster_length; i++) {
        curr_atom = cluster[i].split(/\s+/);

        curr_name = curr_atom[0];
        curr_load = parseFloat(curr_atom[1]);
        curr_x = parseFloat(curr_atom[2]) * scale_const - 2 * min_x;
        curr_y = parseFloat(curr_atom[3]) * scale_const - 2 * min_y;
        curr_z = parseFloat(curr_atom[4]) * scale_const - 2 * min_z;

        var atom = new THREE.Mesh(atom_geometry, new THREE.MeshLambertMaterial({
            color: curr_load / 100 * curr_load / 100 * 0xffffff
        }));

        atom.position.x = curr_x;
        atom.position.y = curr_y;
        atom.position.z = curr_z;
        atom.scale.x = Math.log(curr_load + size_eq) * size_const;
        atom.scale.y = Math.log(curr_load + size_eq) * size_const;
        atom.scale.z = Math.log(curr_load + size_eq) * size_const;
        atom.userData.id = (i + 1).toString();
        atom_id[i] = [atom.userData.id, parseFloat(curr_atom[2]), parseFloat(curr_atom[3]), parseFloat(curr_atom[4]), curr_atom[0]];
        scene.add(atom);
    }

    for (var i = 0; i < cluster_length; i++) {
        for (var j = 0; j < cluster_length; j++) {
            if (bond_array[i][j][0] != 0 && bond_array[i][j][0] < bond_item.Bond_Distance) {
                var bond = new THREE.Mesh(bond_geometry, new THREE.MeshLambertMaterial({
                    color: 0.9 * 0xffffff
                }));
                bond.position.x = bond_array[i][j][4];
                bond.position.y = bond_array[i][j][5];
                bond.position.z = bond_array[i][j][6];
                bond.scale.y = bond_array[i][j][0];
                var vector = new THREE.Vector3(bond_array[i][j][7], bond_array[i][j][8], bond_array[i][j][9]);
                var axis = new THREE.Vector3(0, 1, 0);
                bond.quaternion.setFromUnitVectors(axis, vector.clone().normalize());
                scene.add(bond);

            }
        }
    }
    max_dist_bond = bond_item.Bond_Distance;

}


function make_CP_prop() {

    source_points = new Array((source_raw.length - 5) / 3)
    var count = 1;
    for (var i = 0; i < source_raw.length; i++) {
        
        if (i == 0) {
            var source_atoms_temp = source_raw[i].split(/\s+/);
            source_contrib = new Array((source_atoms_temp.length - 4) / 2);
            for (var j = 0; j < source_atoms_temp.length; j++) {
                if (j > 0 && j < source_atoms_temp.length - 3 && j % 2 != 0) {
                    source_contrib[j - 1] = new Array((source_raw.length - 5) / 3 + 1)
                    source_contrib[j - 1][0] = source_atoms_temp[j] + source_atoms_temp[j + 1];
                }
            }
            source_contrib.clean(undefined);
        } else {
            var temp = source_raw[i].split(/\s+/);
            if (((temp[0] == "ACP") || (temp[0] == "BCP")) && ((temp[2] == "%") || (temp[3] == "%"))) {
                if (temp[2] == "%") {
                    source_points[count] = temp[0] + " " + temp[1];
                    for (var k = 3; k < temp.length - 2; k++) {
                        source_contrib[k - 3][count] = temp[k];

                    }
                } else {
                    source_points[count] = temp[0] + " " + temp[1] + " " + temp[2];
                    for (var k = 4; k < temp.length - 2; k++) {
                        source_contrib[k - 4][count] = temp[k];
                    }
                }
                count++;
            }
        }
    }
    source_points.clean(undefined);

    CPproperties = new Array(CPraw.length);
    var count2 = 0;
    for (var i = 0; i < CPraw.length; i++) {
        if (CPraw[i].split(/\s+/)[1] == "Position") {
            CPproperties[count2] = new Array(11+source_points.length);
            CPproperties[count2][2] = [parseFloat(CPraw[i].split(/\s+/)[3]), parseFloat(CPraw[i].split(/\s+/)[4]), parseFloat(CPraw[i].split(/\s+/)[5])];
            CPproperties[count2][3] = [parseFloat(CPraw[i + 1].split(/\s+/)[5])];
            CPproperties[count2][4] = [parseFloat(CPraw[i + 5].split(/\s+/)[5])];
            CPproperties[count2][5] = [parseFloat(CPraw[i + 7].split(/\s+/)[5])];
            CPproperties[count2][6] = [parseFloat(CPraw[i + 8].split(/\s+/)[5])];
            CPproperties[count2][7] = [parseFloat(CPraw[i + 9].split(/\s+/)[6])];
            CPproperties[count2][8] = [parseFloat(CPraw[i + 10].split(/\s+/)[5])];
            CPproperties[count2][9] = [parseFloat(CPraw[i + 11].split(/\s+/)[5])];
            count2++;
        }
    }
    CPproperties.clean(undefined);


    for (var j = 0; j < CPproperties.length; j++) {
        var atom_test = get_atom_id(CPproperties[j][2][0], CPproperties[j][2][1], CPproperties[j][2][2], 0.25)
        if (atom_test[0] != 0) {
            CPproperties[j][0] = atom_test[0];
            CPproperties[j][1] = atom_test[1];
            for(var k = 0; k < source_contrib.length; k++){
                if(source_contrib[k][0] == atom_test[1]){
                    for(var l = 1; l < source_contrib[k].length; l++){
                        CPproperties[j][9+l] = source_contrib[k][l]
                    }
                }
            }
        }
    }
}


function onDocumentMouseMove(event) {
    event.preventDefault();
    mouse.x = (event.clientX / window.innerWidth) * 2 - 1;
    mouse.y = -(event.clientY / window.innerHeight) * 2 + 1;
}

function get_atom_id(x, y, z, interval) {
    for (var i = 0; i < atom_id.length; i++) {
        if (Math.abs(x - atom_id[i][1]) < interval && Math.abs(y - atom_id[i][2]) < interval && Math.abs(z - atom_id[i][3]) < interval) {
            return [atom_id[i][0], atom_id[i][4]];
        }
    }
    return 0;
}

function make_GUI() {
    cp_nam_item = new cp_nam_item_fun();
    cp_pos_item = new cp_pos_item_fun();
    cp_den_item = new cp_den_item_fun();
    cp_lag_item = new cp_lag_item_fun();
    cp_ham_item = new cp_ham_item_fun();
    cp_pot_item = new cp_pot_item_fun();
    cp_end_item = new cp_end_item_fun();
    cp_lap_item = new cp_lap_item_fun();
    cp_loc_item = new cp_loc_item_fun();
    source_item = new source_item_fun();
    bond_item = new bond_item_fun();
    source_index_item = new source_index_item_fun();
    source_name_item = new source_name_item_fun();
    gui = new dat.GUI();
    gui.add(cp_nam_item, 'Atom_Name').listen();
    gui.add(cp_pos_item, 'Atom_Position_Bohr').listen();
    var f1 = gui.addFolder('Various Properties')
    f1.add(cp_den_item, 'Electron_Density').listen();
    f1.add(cp_lag_item, 'Lagrangian_Kinetic_Energy').listen();
    f1.add(cp_ham_item, 'Hamiltonian_Kinetic_Energy').listen();
    f1.add(cp_pot_item, 'Potential_Energy_Density').listen();
    f1.add(cp_end_item, 'Energy_Density').listen();
    f1.add(cp_lap_item, 'Laplacian_Of_Electron_Density').listen();
    f1.add(cp_loc_item, 'Electron_Localization_Function').listen();
    var f2 = gui.addFolder('Source Contributions')
    f2.add(source_item, 'Source_Contribution').listen();    
    f2.add(source_index_item, 'IP_Index', 1, source_points.length);
    f2.add(source_name_item, 'IP_Name').listen();
    gui.add(bond_item, 'Bond_Distance', 0, 10);
}

function update_GUI(cp_nam, cp_pos, cp_den, cp_lag, cp_ham, cp_pot, cp_end, cp_lap, cp_loc, source, source_name) {

    cp_nam_item.Atom_Name = cp_nam.toString();
    cp_pos_item.Atom_Position_Bohr = parseFloat(cp_pos[0]).toFixed(4).toString() + ", " + parseFloat(cp_pos[1]).toFixed(4).toString() + ", " + parseFloat(cp_pos[2]).toFixed(4).toString();
    cp_den_item.Electron_Density = parseFloat(cp_den).toFixed(4).toString();
    cp_lag_item.Lagrangian_Kinetic_Energy = parseFloat(cp_lag).toFixed(4).toString();
    cp_ham_item.Hamiltonian_Kinetic_Energy = parseFloat(cp_ham).toFixed(4).toString();
    cp_pot_item.Potential_Energy_Density = parseFloat(cp_pot).toFixed(4).toString();
    cp_end_item.Energy_Density = parseFloat(cp_end).toFixed(4).toString();
    cp_lap_item.Laplacian_Of_Electron_Density = parseFloat(cp_lap).toFixed(4).toString();
    cp_loc_item.Electron_Localization_Function = parseFloat(cp_loc).toFixed(4).toString();
    source_item.Source_Contribution = parseFloat(source).toFixed(4).toString();
    source_name_item.IP_Name = source_name;

};

Array.prototype.clean = function (deleteValue) {
    for (var i = 0; i < this.length; i++) {
        if (this[i] == deleteValue) {
            this.splice(i, 1);
            i--;
        }
    }
    return this;
};