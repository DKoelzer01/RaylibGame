// survival/creative-invulnerable (#1561): creative mode takes no damage of any kind and never dies. In creative: a 500-point
// hit, an empty oxygen tank for 3 s (suffocation) and a 40 m fall take nothing and the player lives. Then survival, as a
// control: the same kind of hit does take health.
#include "scenario.h"
#include "combat.h"
#include "player.h"
#include "playersession.h"
#include "survival.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>

namespace {
float health0 = 0.0f;
DamageResult Hit(float amount) {
    DamageEvent e; e.amount = amount; e.kind = DamageKind::Melee; e.hasSource = true;
    e.source = player.position + DVec3(Vector3Scale(player.forward, 2.0f));
    return Combat_Damage(LocalSession(), e);
}
}

SCENARIO("survival/creative-invulnerable") {
    s.tags = { "skills" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;
    s.step("settle, creative on", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        terrainTool.creative = true;
        player.flyMode = false;
        survival.health = 100.0f;
        return true;
    });
    s.step("creative: a 500-point hit and 3 s without oxygen take nothing", 10.0, [](ScenarioCtx& c) {
        if (c.t < 0.1) { health0 = survival.health; Hit(500.0f); }
        survival.oxygen = 0.0f; // (suffocating, every step)
        if (c.t < 3.0) return false;
        survival.oxygen = 100.0f;
        c.check("creative: still alive", !survival.dead);
        c.near("creative: health unchanged", survival.health, health0, 0.01f);
        return true;
    });
    s.step("creative: a 40 m fall takes nothing", 20.0, [](ScenarioCtx& c) {
        static bool dropped = false;
        if (!dropped) {
            const Vector3 up = Vector3Normalize(c.body->toLocal((player.position - c.world.renderOrigin).toFloat()));
            player.flyMode = false;
            player.position = player.position + DVec3(Vector3Scale(c.body->dirToWorld(up), 40.0f));
            player.velocity = Vector3{ 0, 0, 0 };
            health0 = survival.health;
            dropped = true;
            return false;
        }
        c.camera(Vector3Add(c.body->toLocal((player.position - c.world.renderOrigin).toFloat()), { 6.0f, 2.0f, 6.0f }),
                 c.body->toLocal((player.position - c.world.renderOrigin).toFloat()), Vector3{ 0, 1, 0 });
        if (!player.grounded || c.t < 1.0) return false;
        c.check("creative: landed alive", !survival.dead);
        c.near("creative: no fall damage", survival.health, health0, 0.01f);
        return true;
    });
    s.step("survival (control): the same hit takes health", 10.0, [](ScenarioCtx& c) {
        if (c.t < 0.1) { terrainTool.creative = false; return false; } // (the flag follows the mode on the next step)
        if (c.t < 0.3) return false;
        health0 = survival.health;
        Hit(30.0f);
        c.lt("survival: a 30-point hit takes health", survival.health, health0 - 1.0f);
        return true;
    });
}
