//  This file is part of Project Name
//  Copyright (C) Michigan State University, 2017.
//  Released under the MIT Software license; see doc/LICENSE

#include "web/web.h"
#include "web/color_map.h"
#include "../memic_model.h"
#include "config/config_web_interface.h"

namespace UI = emp::web;

class HCAWebInterface : public UI::Animate, public HCAWorld{
  // friend class HCAWorld;
  // friend class UI::Animate;

  using color_fun_t = std::function<std::string(int)>;

  MemicConfig config;
  emp::ConfigWebUI config_ui;
  emp::Random r;
  // HCAWorld world;

  // UI::Animate anim;
  UI::Document oxygen_area;
  UI::Document cell_area;
  UI::Document controls;
  UI::Canvas oxygen_display;
  UI::Canvas cell_display;
  // UI::Canvas clade_display;
  const double display_cell_size = 5;

  UI::Button toggle;
  UI::Style button_style;

  color_fun_t cell_color_fun;
  UI::Selector cell_color_control;
  color_fun_t phylo_depth_color_fun = [this](int cell_id) {
                                        auto taxon = systematics[0].DynamicCast<emp::Systematics<Cell, int>>()->GetTaxonAt(cell_id);
                                        double depth = taxon->GetDepth();
                                        double max_depth = systematics[0]->GetMaxDepth();
                                        double depth_hue = 0;
                                        if (max_depth > 0) {
                                          depth_hue = depth * 280.0/max_depth;
                                        }
                                        return emp::ColorHSL(depth_hue,50,50);
                                     };
  color_fun_t hif1alpha_color_fun = [this](int cell_id) {
                                        double hue = pop[cell_id]->hif1alpha * 280.0;
                                        return emp::ColorHSL(hue,50,50);
                                     };

  public:
  HCAWebInterface() : config_ui(config), oxygen_area("oxygen_area"), cell_area("cell_area"), controls("control_area"), 
    oxygen_display(400, 400, "oxygen_display"), cell_display(400, 400, "cell_display"),
    cell_color_control("cell_color_control")
    // : anim([this](){DoFrame();}, oxygen_display, cell_display) 
  {
    SetupInterface();   
  }

  void SetupInterface() {
    GetRandom().ResetSeed(config.SEED());
    Setup(config, true);
    oxygen_display.SetSize(WORLD_X * display_cell_size, WORLD_Y * display_cell_size);
    oxygen_display.Clear("black");

    cell_display.SetSize(WORLD_X * display_cell_size, WORLD_Y * display_cell_size);
    cell_display.Clear("black");

    oxygen_area << "<h1 class='text-center'>Oxygen</h1>" << oxygen_display;
    cell_area << "<h1 class='text-center'>Cells</h1>" << cell_display;

    cell_color_fun = phylo_depth_color_fun;

    cell_color_control.SetOption("Phylogenetic Depth", 
                                 [this](){
                                     cell_color_fun = phylo_depth_color_fun;
                                     RedrawCells();
                                 }, 0);

    cell_color_control.SetOption("Hif1-alpha stain", 
                                 [this](){
                                     cell_color_fun = hif1alpha_color_fun;
                                     RedrawCells();
                                 }, 1);


    toggle = GetToggleButton("but_toggle");
    button_style.AddClass("btn");
    button_style.AddClass("btn-primary");
    toggle.SetCSS(button_style);
    controls << toggle;
    controls << " " << cell_color_control;

    oxygen_display.On("click", [this](int x, int y){OxygenClick(x, y);});;
    RedrawOxygen();
    RedrawCells();
    
    emp::web::OnDocumentReady([](){
        // EM_ASM(d3.select("#but_toggle").classed("btn btn-primary", true););
        EM_ASM($('select').selectpicker('setStyle', 'btn-primary'););
    });

    config_ui.SetOnChangeFun([this](const std::string & val){ std::cout << "New val: " << val<<std::endl;;InitConfigs(config);});
    config_ui.Setup();
    controls << config_ui.GetDiv();
  }

  void DoFrame() {
    // std::cout << frame_count << " " << GetStepTime() << std::endl;
    UpdateOxygen();

    if (frame_count % DIFFUSION_STEPS_PER_TIME_STEP == 0) {
      RunStep();
      RedrawOxygen();
      RedrawCells();
    }

  }

  void RedrawOxygen() {
    // oxygen_display.SetSize(WORLD_X * display_cell_size, WORLD_Y * display_cell_size);
    oxygen_display.Freeze();
    oxygen_display.Clear("black");

    for (int x = 0; x < WORLD_X; x++) {
      for (int y = 0; y < WORLD_Y; y++) {
        double o2 = oxygen->GetVal(x,y);

        o2 *= 360;
        if (o2 > 360) {
          o2 = 360;
        } else if (o2 < 0) {
          o2 = 0;
        }
        std::string color = emp::ColorHSL(o2,50,50);
        oxygen_display.Rect(x*display_cell_size, y*display_cell_size, display_cell_size, display_cell_size, color, color);
      }
    }
    oxygen_display.Activate();
  }

  void RedrawCells(){
    // cell_display.SetSize(WORLD_X * display_cell_size, WORLD_Y * display_cell_size);
    cell_display.Freeze();
    cell_display.Clear("black");

    for (int x = 0; x < WORLD_X; x++) {
      for (int y = 0; y < WORLD_Y; y++) {
        int cell_id = x + y * WORLD_X;
        if (IsOccupied(cell_id)) {
          // auto taxon = systematics[0].DynamicCast<emp::Systematics<Cell, int>>()->GetTaxonAt(cell_id);
          // double depth = taxon->GetDepth();
          // double depth_hue = depth * 280.0/systematics[0]->GetMaxDepth();
          // std::cout << depth_hue << std::endl;
          std::string color = cell_color_fun(cell_id);

          // double clade_hue = pop[cell_id]->clade;
          // // std::cout << "Coloring: " << clade_hue << " " <<clade_hue/next_clade << std::endl;
          // clade_hue *= 260.0/next_clade;
          // std::string color = emp::ColorHSL(clade_hue,50,50);
          cell_display.Rect(x*display_cell_size, y*display_cell_size, display_cell_size, display_cell_size, color, color);
 
          // switch(GetOrg(cell_id).state) {
          //   case CELL_STATE::HEALTHY:
          //     cell_display.Rect(x*display_cell_size, y*display_cell_size, display_cell_size, display_cell_size, "green","green");
          //     break;
          //   case CELL_STATE::TUMOR:
          //     cell_display.Rect(x*display_cell_size, y*display_cell_size, display_cell_size, display_cell_size, "blue", "blue");
          //     break;
          //   default:
          //     std::cout << "INVALID CELL STATE" << std::endl;
          //     break;
          // }
        }        
      }
    }

    cell_display.Activate();

  }

  void OxygenClick(int x, int y) {
    
    // std::cout << "x: " << in_x << " y: " << in_y  <<std::endl;
    // double x = canvas.GetAdjustedX(in_x);
    // double y = canvas.GetAdjustedY(in_y);

    // UI::Canvas canvas = doc.Canvas("world_canvas");
    const double canvas_x = (double) oxygen_display.GetWidth();
    const double canvas_y = (double) oxygen_display.GetHeight();
    double px = ((double) x) / canvas_x;
    double py = ((double) y) / canvas_y;

    size_t pos_x = (size_t) (WORLD_X * px);
    size_t pos_y = (size_t) (WORLD_Y * py);
    // std::cout << "x: " << x << " y: " << y << "WORLD_X: " << WORLD_X << " WORLD_Y: " << WORLD_Y << " canvas_x: " << canvas_x <<" canvas_y: " << canvas_y  << " px: " << px <<  " py: " << py <<" pos_x: " << pos_x << " pos_y: " << pos_y <<std::endl;
    oxygen->SetNextVal(pos_x, pos_y, 1);
    oxygen->SetVal(pos_x, pos_y, 1);
  }

};


HCAWebInterface interface;

int main()
{
}
