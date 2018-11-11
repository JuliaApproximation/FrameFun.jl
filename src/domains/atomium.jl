# atomium.jl

###
# The atomium: a famous building in Belgium
###

function atomium()
    sphere1 = ball(0.25)
    spheres = union(sphere1,sphere1 + SVector( 0.6, 0.6, 0.6))
    spheres= union(spheres, sphere1 + SVector( 0.6, 0.6,-0.6))
    spheres= union(spheres, sphere1 + SVector( 0.6,-0.6, 0.6))
    spheres= union(spheres, sphere1 + SVector( 0.6,-0.6,-0.6))
    spheres= union(spheres, sphere1 + SVector(-0.6, 0.6, 0.6))
    spheres= union(spheres, sphere1 + SVector(-0.6, 0.6,-0.6))
    spheres= union(spheres, sphere1 + SVector(-0.6,-0.6, 0.6))
    spheres= union(spheres, sphere1 + SVector(-0.6,-0.6,-0.6))
    cyl1 = cylinder(0.10, 1.2)
    spheres= union(spheres, cyl1 + SVector(0.6, 0.6, -0.6));
    spheres= union(spheres, cyl1 + SVector(0.6,-0.6, -0.6));
    spheres= union(spheres, cyl1 + SVector(-0.6, 0.6,-0.6));
    spheres= union(spheres, cyl1 + SVector(-0.6,-0.6,-0.6));
    cyl2 = DomainSets.rotate(cyl1, 0.0, pi/2.0, 0.0)
    spheres= union(spheres, cyl2 + SVector(-0.6,  0.6, 0.6))
    spheres= union(spheres, cyl2 + SVector(-0.6, -0.6, 0.6))
    spheres= union(spheres, cyl2 + SVector(-0.6,  0.6,-0.6))
    spheres= union(spheres, cyl2 + SVector(-0.6, -0.6,-0.6))
    cyl2b = DomainSets.rotate(cyl1, pi/2.0, 0.0, 0.0)
    spheres= union(spheres, cyl2b + SVector( 0.6, 0.6,-0.6))
    spheres= union(spheres, cyl2b + SVector(-0.6, 0.6, 0.6))
    spheres= union(spheres, cyl2b + SVector( 0.6, 0.6, 0.6))
    spheres= union(spheres, cyl2b + SVector(-0.6, 0.6,-0.6))
    cyl3 = cylinder(0.10, 1.2*sqrt(3))
    cyl3 = DomainSets.rotate(cyl3, 0.0, acos(1/sqrt(3)), 0.0)
    cyl3 = DomainSets.rotate(cyl3, 0.0, 0.0, pi/4)
    spheres= union(spheres, cyl3 + SVector( -0.6, -0.6, -0.6))
    cyl4 = cylinder(0.10, 1.2*sqrt(3))
    cyl4 = DomainSets.rotate(cyl4, 0.0, -acos(1/sqrt(3)), 0.0)
    cyl4 = DomainSets.rotate(cyl4, 0.0, 0.0, pi/4)
    spheres= union(spheres, cyl4 + SVector(  0.6, 0.6,-0.6))
    cyl5 = cylinder(0.10, 1.2*sqrt(3))
    cyl5 = DomainSets.rotate(cyl5, 0.0, acos(1/sqrt(3)), 0.0)
    cyl5 = DomainSets.rotate(cyl5, 0.0, 0.0, -pi/4)
    spheres= union(spheres, cyl5 + SVector( -0.6, +0.6, -0.6))
    cyl6 = cylinder(0.10, 1.2*sqrt(3))
    cyl6 = DomainSets.rotate(cyl6, 0.0, -acos(1/sqrt(3)), 0.0)
    cyl6 = DomainSets.rotate(cyl6, 0.0, 0.0, -pi/4)
    spheres= union(spheres, cyl6 + SVector( 0.6, -0.6, -0.6))
    atomium = spheres
end
