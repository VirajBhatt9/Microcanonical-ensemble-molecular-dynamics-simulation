// constants declaration
sig = 1;
ep = 1;
m = 1;
lx = 20;
ly = lx;
lz = lx;
np = 512;
ni = 2000;
dt = 0.005;
va = sqrt(12);
rc = 2.5 * sig;
avx = 0;
avy = 0;
avz =0;
t2 = (dt * dt) / (2 * m);
t1 = dt / (2 * m);
fc = (24 / (rc^7)) * ((2 / (rc^6)) - 1);
Vc = (4 / (rc^6)) * ((1/(rc^6)) -1);

// function to Initialie Positions
function [x, y, z] = posinit()
    p1 = 2; dp =2; p2 =16, kk =1;
    for p = p1: dp: p2
        for q = p1: dp: p2
            for r = p1: dp: p2
                x(kk) = p;
                y(kk) = q;
                z(kk) = r;
                kk = kk + 1;
            end
        end
    end        
endfunction

// function to Initialize Velocities
function [vx, vy, vz] = velinit()
    for i = 1: np
        vx(i) = va*(rand() - 0.5);
        vy(i) = va*(rand() - 0.5);
        vz(i) = va*(rand() - 0.5);    
    end
    for i = 1: np
        avx = vx(i) + avx;
        avy = vy(i) + avy;
        avz = vz(i) + avz;
    end
    avx = avx / np;
    avy = avy / np;
    avz = avz / np;
    for i = 1: np
        vx(i) = vx(i) - avx;
        vy(i) = vy(i) - avy;
        vz(i) = vz(i) - avz;
    end
endfunction

// function to initialize Forces
function [fx, fy, fz, PEN] = forcecalc(x, y, z)
    fx = zeros(np,1);
    fy = zeros(np,1);
    fz = zeros(np,1);
    PEN = 0;
    for i =1: np - 1
        x1 = x(i);
        y1 = y(i);
        z1 = z(i);
        for j = i + 1: np
            x2 = x(j);
            y2 = y(j);
            z2 = z(j);
            dx = x1 - x2;
            dy = y1 - y2;
            dz = z1 - z2;
            if abs(dx) > (lx / 2)
                dx = (lx - abs(dx)) * (-dx) / abs(dx);
            end 
            if abs(dy) > (ly / 2)
                dy = (ly - abs(dy)) * (-dy) / abs(dy);
            end
            if abs(dz) > (lz / 2)
                dz = (lz - abs(dz)) * (-dz) / abs(dz);
            end
            r = sqrt((dx * dx) + (dy * dy) + (dz * dz))
            if r <= rc
                f = (24 / (r^7)) * ((2 / (r^6)) - 1);
                fm = f - fc;
                fx(i) = fx(i) + fm * (dx / r);
                fy(i) = fy(i) + fm * (dy / r);
                fz(i) = fz(i) + fm * (dz / r);
                fx(j) = fx(j) - fm * (dx / r);
                fy(j) = fy(j) - fm * (dy / r);
                fz(j) = fz(j) - fm * (dz / r);
                V = (4 / (r^6)) * ((1/(r^6)) -1);
                Vm = V - Vc + (r * fc) - (rc * fc);
                PEN = PEN + Vm;
            end
        end
    end
endfunction

// function to Update Positions
function [x, y, z] = updatepos(x,y,z, vx,vy,vz, fx,fy,fz)
    for i = 1: np
        x(i) = x(i) + (vx(i) * dt) + (fx(i) * t2);
        y(i) = y(i) + (vy(i) * dt) + (fy(i) * t2);
        z(i) = z(i) + (vz(i) * dt) + (fz(i) * t2);
        if x(i) > lx
            then x(i) = x(i) - lx;
        elseif x(i) < 0
             then x(i) = x(i) + lx;
        end
        if y(i) > ly
            then y(i) = y(i) - ly;
        elseif y(i) < 0
             then y(i) = y(i) + ly;
        end
        if z(i) > lz
            then z(i) = z(i) - lz;
        elseif z(i) < 0
             then z(i) = z(i) + lz;
        end 
    end
endfunction

// function to Update Velocities
function [vx, vy, vz, KEN] = updatevel(vx,vy,vz, fx,fy,fz, fox,foy,foz)
    for i = 1: np
        vx(i) = vx(i) + (fox(i) + fx(i)) * t1;
        vy(i) = vy(i) + (foy(i) + fy(i)) * t1;
        vz(i) = vz(i) + (foz(i) + fz(i)) * t1;
    end
    KEN = 0;
    for i = 1: np
        KEN = KEN + 0.5 * ((vx(i)^2) + (vy(i)^2) + (vz(i)^2))
    end
endfunction

// Main program
[x, y, z] = posinit();
x1 = x;
y1 =y;
z1 = z;
[vx, vy, vz] = velinit();
[fx, fy, fz, PEN] = forcecalc(x, y, z);
for k = 1: ni
    [x, y, z] = updatepos(x,y,z, vx,vy,vz, fx,fy,fz)
    fox = fx;
    foy = fy;
    foz = fz;
    [fx, fy, fz, PEN] = forcecalc(x, y, z);
    [vx, vy, vz, KEN] = updatevel(vx,vy,vz, fx,fy,fz, fox,foy,foz)
    PE(k) = PEN / np;
    KE(k) = KEN / np;
    nit(k) = k;
    disp('Iterations: ' + string(k));
end

subplot(2,2,1)
scatter3d(x1,y1,z1,"fill");
replot([0,0,20,20]);
title('Initial position of particles','fontsize',6);
subplot(2,2,2)
scatter3d(x,y,z,"fill");
title('Final position of particles','fontsize',6);

TE = KE + PE;
subplot(2,2,3)
plot(nit',KE',nit',PE',nit',TE','linewidth',6);
legend('Kinetic Energy', 'Potential Energy', 'Total Energy', 5)
xlabel('No. of interations', 'fontsize', 6);
ylabel('Energy', 'fontsize', 6);

for i = 1:np
    velocity(i) = sqrt((vx(i)^2) + (vy(i)^2) + (vz(i)^2))
end
subplot(2,2,4)
histplot(20, velocity, polygon=%t);
xlabel('Speed of the particles','fontsize',6);
