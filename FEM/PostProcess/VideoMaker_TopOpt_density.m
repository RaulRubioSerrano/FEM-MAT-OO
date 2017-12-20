classdef  VideoMaker_TopOpt_density < VideoMaker_TopOpt
  properties
      
      
      
  end
  
  methods (Access = public)

        function Make_video_design_variable(obj,output_video_name)
            post = Postprocess_TopOpt_density();
            field2print = post.density_name;
            componentfield = post.density_name_component;
            obj.Print_make_video_density_field(field2print,componentfield,output_video_name)
        end
        
        
        function Make_video_design_variable_reg(obj,output_video_name)
            post = Postprocess_TopOpt_density();
            field2print = post.density_name_reg;
            componentfield = post.density_name_component_reg;
            obj.Print_make_video_density_field(field2print,componentfield,output_video_name)
        end
        
        
  end
        
  methods (Access = private)
        function Print_make_video_density_field(obj,field2print,componentfield,output_video_name)
            file_tcl_name = 'tcl_gid.tcl';
            file_list = obj.create_file_list(obj.iterations_to_print,obj.file_name,obj.files_folder);
            file_tcl_name_with_path = fullfile(obj.files_folder,file_tcl_name);
            fid = fopen(file_tcl_name_with_path,'w+');
            fprintf(fid,'GiD_Process PostProcess \n');
            fprintf(fid,['set arg1 "',file_list,'"\n']);
            fprintf(fid,['set arg2 "',output_video_name,'"\n']);
            fprintf(fid,['set arg3 "',field2print,'"\n']);
            fprintf(fid,['set arg4 "',componentfield,'"\n']);
            fprintf(fid,['source "',fullfile(pwd,'FEM','PostProcess','Make_Video_density.tcl'),'"\n']);
            fprintf(fid,['Make_Video_density $arg1 $arg2 $arg3 $arg4 \n']);
            fprintf(fid,['GiD_Process Mescape Quit']);
            fclose(fid);
            obj.execute_tcl_files(obj.gidPath,file_tcl_name_with_path)
        end
        
    end
end