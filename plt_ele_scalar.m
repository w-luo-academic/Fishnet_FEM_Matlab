function [] = plt_ele_scalar(v_his,num_e,connect,coord,m,n,n_frame)
    sig_max = 20; %max(abs(v_his(:,n_frame)));
    for ii = 1:num_e
       node_1 = connect(ii,1);
       node_2 = connect(ii,2);
       plot([coord(node_1,1), coord(node_2,1)],[coord(node_1,2), coord(node_2,2)],...
           'linewidth',2,'Color',[0 0 0] + 1 - abs(v_his(ii,n_frame))/sig_max...
       )
       pbaspect([n m 1])
       set(gca,'XTick',[],'YTick',[])
       hold on;
    end
end