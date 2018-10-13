classdef pop
    % population of Izhikevich neurons
    properties
        Ne
        Ni
        re
        ri
        a
        b
        c
        d
        w
        wmax
        v
        v_md
        v_mp
        v_rekt_md
        v_rekt_mp
        vy
        x
        x0
        u 
        STD
    end
    methods
        function obj=pop(a,b)
            obj.Ne = a;
            obj.Ni = b;
            obj.re=rand(obj.Ne,1);          obj.ri=rand(obj.Ni,1);
            obj.a=[0.02*ones(obj.Ne,1);     0.02+0.08*obj.ri];
            obj.b=[0.2*ones(obj.Ne,1);      0.25-0.05*obj.ri];
            obj.c=[-65+15*obj.re.^2;        -65*ones(obj.Ni,1)];
            obj.d=[8-6*obj.re.^2;           2*ones(obj.Ni,1)];
            obj.w=[normrnd(2.5,.75,obj.Ne+obj.Ni,obj.Ne),    normrnd(-2.5,.75,obj.Ne+obj.Ni,obj.Ni)];
            obj.wmax = 10;
            obj.v=-65*ones(obj.Ne+obj.Ni,1);          % Initial values of v
            obj.v_md = obj.v;               % low pass of v for the depression term
            obj.v_mp = obj.v;               % low pass of v for the potentiation term
            obj.v_rekt_md = zeros(obj.Ne+obj.Ni,1);      % low-pass of v thresholded for the depression term
            obj.v_rekt_mp = zeros(obj.Ne+obj.Ni,1);      % low-pass of v thresholded for the potentiation term
            obj.vy = zeros(obj.Ne+obj.Ni,1);             % voltage thresholded
            obj.x = zeros(obj.Ne+obj.Ni,obj.Ne);             % presyn low pass
            obj.x0 = obj.x;
            obj.u = obj.b.*obj.v; 
            obj.STD = ones(obj.Ne+obj.Ni,obj.Ne);
        end
    end
end