% Simulación Avanzada de Dinámica Molecular en 2D con Termodinámica
clear; clc; close all;

% % % % % % --- 1. PARÁMETROS DE LA SIMULACIÓN ---
N = 36;                % Número de partículas (cuadrado perfecto)
L = 12;              % Longitud de la caja de simulación
dt = 0.001;            % Paso de tiempo (time step)
steps = 1500;           % Número total de pasos
sigma = 1.0;           % Parámetro sigma de Lennard-Jones
epsilon = 1.0;         % Parámetro epsilon de Lennard-Jones
mass = 1.0;            % Masa de las partículas
cutoff = 3.0 * sigma;  % Radio de corte

% --- 2. INICIALIZACIÓN ---
n_side = round(sqrt(N));
spacing = L / (n_side + 1);
[x, y] = meshgrid(spacing:spacing:L-spacing, spacing:spacing:L-spacing);
pos = [x(:), y(:)];
pos = pos(1:N, :);

% Velocidades iniciales (temperatura inicial)
vel = randn(N, 2) * 0.8;
vel = vel - mean(vel, 1); % Eliminar velocidad del centro de masa

acc = zeros(N, 2);

% Vectores para almacenar series temporales (para graficar termodinámica)
time_vec = zeros(steps, 1);
temp_vec = zeros(steps, 1);
energy_vec = zeros(steps, 1);

% --- 3. BUCLE PRINCIPAL ---
figure('Color', 'white', 'Position', [50, 50, 1100, 550]);

for step = 1:steps
    % Velocity-Verlet: Paso 1 (Posiciones)
    pos = pos + vel * dt + 0.5 * acc * dt^2;
    pos = mod(pos, L); % Condiciones de Frontera Periódicas
    
    acc_old = acc;
    acc = zeros(N, 2);
    potential_energy = 0;
    
    % Cálculo de fuerzas y energía potencial
    for i = 1:N
        for j = i+1:N
            dx = pos(i,1) - pos(j,1);
            dy = pos(i,2) - pos(j,2);
            
            % Imagen mínima
            dx = dx - L * round(dx / L);
            dy = dy - L * round(dy / L);
            
            r2 = dx^2 + dy^2;
            if r2 < cutoff^2
                r6 = (sigma^2 / r2)^3;
                r12 = r6^2;
                
                potential_energy = potential_energy + 4 * epsilon * (r12 - r6);
                f_scalar = 24 * epsilon * (2 * r12 - r6) / r2;
                
                acc(i, :) = acc(i, :) + (f_scalar * [dx, dy]) / mass;
                acc(j, :) = acc(j, :) - (f_scalar * [dx, dy]) / mass;
            end
        end
    end
    
    % Velocity-Verlet: Paso 2 (Velocidades)
    vel = vel + 0.5 * (acc_old + acc) * dt;
    
    % --- CÁLCULO DE PROPIEDADES TERMODINÁMICAS ---
    kinetic_energy = 0.5 * mass * sum(sum(vel.^2));
    total_energy = kinetic_energy + potential_energy;
    
    % Temperatura instantánea en 2D (con 2N - 2 grados de libertad por el centro de masa)
    temperature = (2 * kinetic_energy) / (2 * N - 2);
    
    % Guardar datos en arreglos
    time_vec(step) = step * dt;
    temp_vec(step) = temperature;
    energy_vec(step) = total_energy
    
    % --- 4. VISUALIZACIÓN DINÁMICA (Doble Panel) ---
    if mod(step, 4) == 0
        % Panel izquierdo: Trayectorias y posiciones
        subplot(1, 2, 1);
        scatter(pos(:,1), pos(:,2), 70, 'filled', 'MarkerFaceColor', '#D95319', 'MarkerEdgeColor', 'k');
        xlim([0 L]); ylim([0 L]);
        title(sprintf('Dinámica Molecular 2D - Paso: %d', step), 'FontSize', 11);
        xlabel('Posición X'); ylabel('Posición Y');
        grid on; box on;
        
        % Panel derecho: Evolución termodinámica (Temperatura y Energía)
        subplot(1, 2, 2);
        plot(time_vec(1:step), temp_vec(1:step), 'r-', 'LineWidth', 1.5);
        hold on;
        plot(time_vec(1:step), energy_vec(1:step) / N, 'b--', 'LineWidth', 1.2);
        hold off;
        legend({'Temperatura ($T$)', 'Energía Total / $N$'}, 'Interpreter', 'latex', 'Location', 'best');
        xlim([0 steps*dt]);
        title('Evolución Termodinámica del Sistema', 'FontSize', 11);
        xlabel('Tiempo Simulado (t)'); ylabel('Magnitudes Reducidas');
        grid on; box on;
        
        drawnow;
    end
end