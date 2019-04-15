//  This file is part of Project Name
//  Copyright (C) Michigan State University, 2017.
//  Released under the MIT Software license; see doc/LICENSE

#include "web/web.h"
#include "web/color_map.h"
#include "../memic_model.h"

namespace UI = emp::web;

class HCAWebInterface : public UI::Animate, public HCAWorld{
  // friend class HCAWorld;
  // friend class UI::Animate;

  MemicConfig config;
  emp::Random r;
  // HCAWorld world;

  // UI::Animate anim;
  UI::Document doc;
  UI::Div display_area;
  UI::Canvas oxygen_display;
  UI::Canvas cell_display;
  // UI::Canvas clade_display;
  const double display_cell_size = 5;

  public:
  HCAWebInterface() : doc("emp_base"), display_area("display_area"), oxygen_display(400, 400, "oxygen_display"), cell_display(400, 400, "cell_display")
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

    display_area << oxygen_display << " ";
    display_area << cell_display;
    doc << display_area;
    doc << GetToggleButton("but_toggle");

    oxygen_display.On("click", [this](int x, int y){OxygenClick(x, y);});;
    RedrawOxygen();
    RedrawCells();
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
          double clade_hue = pop[cell_id]->clade;
          // std::cout << "Coloring: " << clade_hue << " " <<clade_hue/next_clade << std::endl;
          clade_hue *= 260.0/next_clade;
          std::string color = emp::ColorHSL(clade_hue,50,50);
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
