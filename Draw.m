clear;
clc;

% Считываем данные из файла
data = csvread('Results.txt');

N = data(1);
M = data(2);
x_left = data(3);
x_right = data(4);
T = data(5);

h = (x_right - x_left)/N;   % Вычисляем величину шага по координате
x = x_left:h:x_right;       % Формируем сетку по координате
tau = (T - 0)/M;            % Вычисляем величину шага по времени
t = 0:tau:T;                % Формируем сетку по времени

counter = 5; % Сдвиг в массиве data, относительно которого начинаем считывать решение
u = zeros(M+1,N+1);
for m = 1:M+1
    for n = 1:N+1
        counter = counter + 1;
        u(m,n) = data(counter);
    end
end

% Рисуем решение
k = 4; % В видео записывается каждый k-ый кадр


    
for m = 1:k:(M+1)
        
    % Рисуем начальное условие
    plot(x,u(1,:),'-og','MarkerSize',3);
    hold on;
    % Рисуем сетку по x
    plot(x,x-x,'-ok','MarkerSize',3);
    hold on;
    %m = 300;
    % Рисуем решениев момент времени t(m + 1)
    plot(x,u(m,:),'-or','MarkerSize',2,'LineWidth',1);
    text(0.05,1.95,['t_{' num2str(m-1,'%4.0f'),'} = ',num2str(t(m),'%4.2f')],'Rotation',0,'Color','k');
    hold on;
    
    hold off;
    axis([x_left x_right (u(1,1) - 0.5) (u(M+1,N+1) + 0.5)]);
    xlabel('x \in (0,1)');
    ylabel('u');
    title('Solution of the Burgers Equation');
    
    drawnow;    
    pause(0.1);
    
    %m
    
    % Сохраняем динамическую картинку в файл, 
    % который может быть преобразован в видео
    %mov((m-1)/k+1) = getframe;
    
end

% Сохраняем видеофайл с презентацией результатов вычислений
%movie2avi(mov, 'BurgersEquation.avi', 'compression', 'None');