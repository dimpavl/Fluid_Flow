function x = granitsy(i,j,dx,L,H, image)
      x = 0;
      %����� 1, ������ 0
      %��������.
      x = image(i,j);
      return;
%������ �������
%    if (i*dx<=L && i*dx>=0 && j*dx<=H && j*dx>=0 && i~=0 && j~=0)
%        x = 1;
%        return;
%    else
%        x = 0;
%        return;
%    end;


   
    
% %   ������� 0.2*H  0.5*L �� ������� ����� � �����
% %   ---|
% %      |--
% %    
% %      |--
% %   ---|
%     if (j*dx<=H*0.4 && j*dx>=0 && j~=0 && i*dx<=L && i*dx>=0 && i~=0)
%         x = 1;
%         return;
%     end;
%     if (j*dx>=H*0.4 && j*dx<=H && j~=0)
%         if (i*dx>=L/4 && i*dx<=3*L/4)
%             x = 1 ;
%             return;
%         else 
%             x = 0;
%             return;            
%     end;
   
%    %������� 0.2*H  0.5*L �� ����� ����� � �������
% %     |----
% %   --|
% %    
% %   --|
% %     |----
%     if (j*dx<=H*0.4 && j*dx>=0 && j~=0 && i*dx>=L/4 && i*dx<=3*L/4)
%         x = 1;
%         return;
%     end;
%     if (j*dx>=H*0.4 && j*dx<=H && j~=0)
%         if (i*dx<=L && i*dx>=0 && i~=0)
%             x = 1 ;
%             return;
%         else 
%             x = 0;
%             return;            
%     end;    
    
%     %������ ������ �������
%     if (j*dx<=H && j*dx>=0 && j~=0)
%         if (i*dx>=L/5*(1-cos(2*pi*j*dx)) && i*dx<=L/5*(4+cos(2*pi*dx*j)))
%             x = 1;
%             return;
%         else
%             x = 0;
%             return;
%         end;
%     end;

% 	 %������ ������ ������� 2
%     if (j*dx<=H && j*dx>=0 && j~=0)
%         if (i*dx>=L/5*(1-cos(2*pi*j*dx)) && i*dx<=L)
%             x = 1;
%             return;
%         else
%             x = 0;
%             return;
%         end;
%     end;
end