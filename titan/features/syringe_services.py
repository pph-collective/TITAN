from . import base_feature
from .. import model as hiv_model
from .. import utils


class SyringeServices(base_feature.BaseFeature):

    name = "syringe_services"

    enrolled_risk = 0.0

    def __init__(self, agent):
        super().__init__(agent)

        self.active = False

    @classmethod
    def update_pop(cls, model: "hiv_model.HIVModel"):
        """
        Update the feature for the entire population (class method).  This is useful for initializing class level trackers that need to be reset each time step, or if enabling a feature for agents needs to be evaluated within the context of the full population (limited slots, or similar).

        Enroll PWID agents in syringe services according to the syring_services timeline and params.

        args:
            model: the instance of HIVModel currently being run
        """
        print(("\n\n!!!!Engaging syringe services program"))
        ssp_num_slots = 0
        ssp_agents = {
            agent for agent in model.pop.pwid_agents if agent.syringe_services.active  # type: ignore[attr-defined]
        }

        for item in model.params.syringe_services.timeline.values():
            if item.start_time <= model.time < item.stop_time:
                cls.enrolled_risk = item.risk

                # linearly interpolate slots between start and stop
                ssp_num_slots = (item.num_slots_stop - item.num_slots_start) / (
                    item.stop_time - item.start_time
                ) * (model.time - item.start_time) + item.num_slots_start

                # If cap indicates all or no agents, do not change
                # otherwise, find true number of slots through distribution
                num_pwid_agents = model.pop.pwid_agents.num_members()
                if 0 < ssp_num_slots < num_pwid_agents:
                    ssp_num_slots = round(
                        model.run_random.betavariate(
                            ssp_num_slots,
                            num_pwid_agents - ssp_num_slots,
                        )
                        * num_pwid_agents
                    )
                break

        target_set = utils.safe_shuffle(
            (model.pop.pwid_agents.members - ssp_agents), model.run_random
        )

        # unenroll agents if above cap
        for agent in ssp_agents.copy():
            if len(ssp_agents) > ssp_num_slots:
                agent.syringe_services.active = False  # type: ignore[attr-defined]
                ssp_agents.remove(agent)
            else:
                break

        # enroll agents if below cap
        if target_set is not None:
            for agent in target_set:
                if len(ssp_agents) < ssp_num_slots:
                    agent.syringe_services.active = True  # type: ignore[attr-defined]
                    ssp_agents.add(agent)
                else:
                    break

        print(
            f"SSP has {ssp_num_slots} target slots with "
            f"{len(ssp_agents)} slots filled"
        )
