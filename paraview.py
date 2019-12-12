RenderView1 = GetRenderView() 

# Turn off "Head Light" 
#RenderView1.LightSwitch = 0 

# Turn off "Light Kit" 
RenderView1.UseLight = 0 

r = GetDisplayProperties() 
r.Ambient = 1.0 
r.Diffuse = 0.0
