function return_data = forward_euler(initial_conditions,time_step,s1,s2,r1,r2)

    %for i in np.arange(steps-1)+1:
    %t[i] = t[i-1] + dt
    %y[i,:] = forwards_euler(e,s,r,y[i-1,:],t,dt,derivative) % 1

    return_data = initial_conditions+time_step+s1+s2+r1+r2;

end