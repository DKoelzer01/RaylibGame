#pragma once

// The Statistics tab (#1199, the factory doc's section 3): items made and used per minute over a window (the last minute,
// hour or 10 hours) with a sparkline, by hand in and out; power made, used and asked for, battery charge; machines by kind
// and status now (the share of time per status on hover); and the time spent in each tier (the tier clock). World scope for
// now (one plot and a team's figures come with the summary sim and the server snapshot).

struct FactoryStatsView;

void StatsTab_Draw();                       // inside the game window's screen (ImGui)
void StatsTab_Refresh();                    // take the view now (it refreshes itself twice a second while shown)
const FactoryStatsView& StatsTab_View();    // (tests) the last view taken
