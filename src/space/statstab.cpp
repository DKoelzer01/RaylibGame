#include "statstab.h"
#include "factory.h"
#include "inventory.h"
#include "plots.h"
#include "survival.h"
#include "imgui.h"

#include <raylib.h>
#include <algorithm>
#include <string>
#include <vector>

namespace {
StatWindow window = StatWindow::Minute;
FactoryStatsView view;
double refreshedAt = -1.0;
bool viewFresh = false;
constexpr double kRefreshSeconds = 0.5; // the view is taken between factory ticks: twice a second is plenty

const char* WindowName(StatWindow w) { return w == StatWindow::Minute ? "last minute" : w == StatWindow::Hour ? "last hour" : "last 10 hours"; }
std::string ItemName(const std::string& id) { const ItemDef* d = itemDB.find(id); return d ? d->name : id; }
std::string Hours(double s) { const int m = (int)(s / 60.0); return TextFormat("%dh %02dm", m / 60, m % 60); }
}

void StatsTab_Refresh() {
    view = Factory::StatsView(window);
    refreshedAt = GetTime();
    viewFresh = true;
}

const FactoryStatsView& StatsTab_View() { return view; }

void StatsTab_Draw() {
    if (!viewFresh || GetTime() - refreshedAt > kRefreshSeconds) StatsTab_Refresh();
    const float s = ImGui::GetFontSize() / 13.0f;
    // The window choice.
    ImGui::TextUnformatted("Window:");
    for (StatWindow w : { StatWindow::Minute, StatWindow::Hour, StatWindow::TenHours }) {
        ImGui::SameLine();
        if (ImGui::RadioButton(WindowName(w), window == w)) { window = w; StatsTab_Refresh(); }
    }
    ImGui::SameLine(0.0f, 20.0f * s);
    if (view.window != window) ImGui::TextDisabled("(not filled yet: showing the %s)", WindowName(view.window));
    else ImGui::TextDisabled("(%d plots, %.1f min of data)", view.plots, view.coveredMinutes);
    ImGui::Separator();

    ImGui::Columns(2, "stats", true);
    // Items: made and used per minute, by hand, and the window's sparkline.
    ImGui::TextUnformatted("Items per minute (hand: put in / taken out)");
    if (ImGui::BeginTable("items", 5, ImGuiTableFlags_RowBg | ImGuiTableFlags_ScrollY | ImGuiTableFlags_SizingStretchProp, ImVec2(0, ImGui::GetContentRegionAvail().y * 0.98f))) {
        ImGui::TableSetupScrollFreeze(0, 1);
        ImGui::TableSetupColumn("Item", ImGuiTableColumnFlags_WidthStretch, 1.8f);
        ImGui::TableSetupColumn("Made", ImGuiTableColumnFlags_WidthStretch, 1.0f);
        ImGui::TableSetupColumn("Used", ImGuiTableColumnFlags_WidthStretch, 1.0f);
        ImGui::TableSetupColumn("Hand", ImGuiTableColumnFlags_WidthStretch, 1.4f);
        ImGui::TableSetupColumn("Trend", ImGuiTableColumnFlags_WidthStretch, 1.6f); // (made over the window)
        ImGui::TableHeadersRow();
        std::vector<const StatLine*> rows;
        for (const StatLine& l : view.items) rows.push_back(&l);
        std::stable_sort(rows.begin(), rows.end(), [](const StatLine* a, const StatLine* b) { return a->madePerMin + a->usedPerMin > b->madePerMin + b->usedPerMin; });
        for (const StatLine* l : rows) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::TextUnformatted(ItemName(l->item).c_str());
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s\nmade %lld, used %lld since the load\nput in %lld, taken out %lld by hand", ItemName(l->item).c_str(), l->made, l->used, l->handIn, l->handOut);
            ImGui::TableNextColumn(); ImGui::Text("%.1f", l->madePerMin);
            ImGui::TableNextColumn(); ImGui::Text("%.1f", l->usedPerMin);
            ImGui::TableNextColumn(); ImGui::Text("%.0f/%.0f", l->handInPerMin, l->handOutPerMin);
            ImGui::TableNextColumn();
            if (!l->madeSpark.empty()) {
                ImGui::PushID(l);
                ImGui::PlotLines("##spark", l->madeSpark.data(), (int)l->madeSpark.size(), 0, nullptr, 0.0f, FLT_MAX, ImVec2(-1, ImGui::GetTextLineHeight()));
                ImGui::PopID();
            }
        }
        if (rows.empty()) { ImGui::TableNextRow(); ImGui::TableNextColumn(); ImGui::TextDisabled("Nothing made or used yet."); }
        ImGui::EndTable();
    }
    ImGui::NextColumn();

    // Power.
    ImGui::TextUnformatted("Power (EU per tick)");
    const double sat = view.euDemandPerTick > 0.0 ? std::min(1.0, view.euUsedPerTick / view.euDemandPerTick) : 1.0;
    ImGui::TextWrapped("Generated %.1f, used %.1f, asked for %.1f", view.euMadePerTick, view.euUsedPerTick, view.euDemandPerTick);
    ImGui::Text("Satisfaction %.0f%%", sat * 100.0);
    if (view.batteryCapacity > 0) ImGui::Text("Batteries %.0f%% (%lld / %lld EU)", 100.0 * view.batteryEnergy / (double)view.batteryCapacity, view.batteryEnergy, view.batteryCapacity);
    for (const StatKind& k : view.kinds)
        if (k.euMadePerTick > 0.0 || k.euUsedPerTick > 0.0)
            ImGui::BulletText("%s: %s %.1f", PartName((PartType)k.type), k.euMadePerTick > 0.0 ? "makes" : "uses", k.euMadePerTick > 0.0 ? k.euMadePerTick : k.euUsedPerTick);
    ImGui::Separator();

    // Machines by kind and status (the bottleneck finder): now, and the share of time since the load.
    ImGui::TextUnformatted("Machines");
    for (const StatKind& k : view.kinds) {
        if (k.count <= 0 || (PartType)k.type == PartType::Chest) continue;
        std::string now;
        for (int st = 0; st < 8; ++st)
            if (k.now[st] > 0) now += TextFormat("%s%d %s", now.empty() ? "" : ", ", k.now[st], MachineStatusName((MachineStatus)st));
        ImGui::Bullet(); ImGui::TextWrapped("%d %s: %s", k.count, PartName((PartType)k.type), now.c_str());
        if (ImGui::IsItemHovered()) {
            long long all = 0;
            for (long long t : k.ticks) all += t;
            std::string share;
            for (int st = 0; st < 8; ++st)
                if (k.ticks[st] > 0) share += TextFormat("%s %.0f%%\n", MachineStatusName((MachineStatus)st), 100.0 * k.ticks[st] / (double)std::max(1LL, all));
            ImGui::SetTooltip("Time since the load:\n%s", share.c_str());
        }
    }
    ImGui::Separator();

    // The tier clock.
    ImGui::TextUnformatted("Time per tier");
    for (int t = 0; t < QuestBook_TierCount(); ++t) {
        const double secs = QuestBook_TierSeconds(t);
        if (secs <= 0.0 && t != QuestBook_Tier()) continue;
        ImGui::BulletText("%s%s: %s", QuestBook_TierName(t), t == QuestBook_Tier() ? " (now)" : "", Hours(secs).c_str());
    }
    ImGui::Columns(1);
}
