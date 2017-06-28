!PATH = expand_path('+' + '/home/shd/IDL_Routines') $
        + path_sep(/search_path) $
        + expand_path('+' + '/home/starfire/STARsoft/IDL') $
        + path_sep(/search_path)
 
; Set up colors for plotting
;device, true_color = 24, retain = 2, decompose = 0
;red = [0,1,1,0,0,1]
;green = [0,1,0,1,0,1]
;blue = [0,1,0,0,1,0]
;if not strcmp(getenv('DISPLAY'),'') then $
;   tvlct, 255*red, 255*green, 255*blue
