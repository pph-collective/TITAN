from typing import List, Dict
from . import base_feature
from .. import agent
from .. import model
from ..parse_params import ObjMap


class PartnerTracing(base_feature.BaseFeature):
    name: str = "partner_tracing"
    stats: List[str] = ["step_1", "step_3", "step_4", "ps_dx", "ps_negative", "ps_prev_dx", "ps_participants", "ps_negative_count", "ps_positive_count", "re_haart", "prev_dx_not_on_haart"]
    
    def __init__(self, agent: "agent.Agent"):
        super().__init__(agent)
        
        self.active = False
        self.time = None
        self.participate = False #willing to participate in PS (i.e name your partners and have them contacted)?
        self.participate_time = None
        self.ps_dx = False #diagnosed with HIV through PS?
        self.ps_dx_time = None #time this agent was diagnosed through PS
        self.tested_negative = False #tested negative through PS?
        self.tested_negative_time = None #time this agent tested negative through PS
        self.info = False # will the agent provide contact info about a partner?
        self.contact = False # will the partner be contacted by PS?
        self.info_stat = 0 #for collecting stats only
        self.contact_stat = 0 #for collecting stats only
        self.ps_prev_dx = False #was this traced agent already previously diagnosed with hiv? (for collecting stats only)
        self.ps_prev_dx_time = None #when did we become aware of the fact that the agent was previously diagnosed? (for collecting stats only)
        self.ps_participant = False #did this partner come back for testing after being contacted? (for collecting stats only)
        self.ps_participant_time = None #when did this agent come back for testing? (for collecting stats only)
        self.re_initiated = False #renitiated in haart after stopping treatment?
        self.prev_dx_not_on_haart = False #is this agent previously diagnosed and is not on haart?
        
    # these stats are used for collecting the number of people who test positive and negative through PS. Note these stats are strictly increasing and should never go down.
    @classmethod
    def init_class(cls, params: "ObjMap"):
        cls.tested_negative_count = 0 
        cls.tested_positive_count = 0
        
    @classmethod    
    def add_negative_count(cls):
        cls.tested_negative_count += 1
    
    @classmethod
    def get_negative_count(cls):
        return cls.tested_negative_count
    
    @classmethod    
    def add_positive_count(cls):
        cls.tested_positive_count += 1
    
    @classmethod
    def get_positive_count(cls):
        return cls.tested_positive_count
    
    def update_agent(self, model: "model.TITAN"):
        """
        Update the agent for this feature for a time step.  Called once per time step in `TITAN.update_all_agents`. Agent level updates are done after population level updates.   Called on only features that are enabled per the params.

        If the agent is was diagnosed last time step, trace their partners. If the agent is traced but not diagnosed, stochastically diagnose.  If the agents tracking has expired, mark them as inactive.

        args:
            model: the instance of TITAN currently being run
        """
        
        params = self.agent.location.params.partner_tracing
        
        if model.time < params.start_time or model.time > params.stop_time:
            return

        agent_exposure = getattr(self.agent, params.exposure)
        
        # if the agent was diagnosed with hiv in the previous time step and is willing to participate, trace their partners
        if (
            agent_exposure.dx
            and agent_exposure.dx_time == model.time - 1
        ):  
            # step 1 
            self.participate = True if model.run_random.random() < params.participate_prob[self.agent.race] else False

            if self.participate:
                self.participate_time = model.time
                
                self.trace_partners(params, model)
                
                        
        # if this agent is traced (i.e. has been contacted by PS) and comes back for testing
        if (
            self.active
            and self.time < model.time
            and model.run_random.random() < params.partner_participate_prob[self.agent.race]
        ):  
            # mark agents who participate
            self.ps_participant = True
            self.ps_participant_time = model.time
            
            # if this traced agent was already previously diagnosed, mark as previously diagnosed (every agent is allowed to enter this if statemnet no more than once)
            if (
                agent_exposure.dx
                and agent_exposure.dx_time < model.time - 1
                and not self.ps_prev_dx # and not already marked
            ):
                self.ps_prev_dx = True
                self.ps_prev_dx_time = model.time
                # is this agent on haart? if not, mark
                if not self.agent.haart.active:
                    self.prev_dx_not_on_haart = True
                    
            # if this traced agent has not been diagnosed, test agent
            if not agent_exposure.dx:
                self.test_agent(agent_exposure, params, model)
                    
                    
            # if agent was traced via PS and was previously diagnosed with hiv and is not on haart, stochastically enroll in haart
            if (
                  self.active 
                  and self.time < model.time
                  and agent_exposure.dx
                  and agent_exposure.dx_time < model.time
                  and not self.agent.haart.active
                  and model.run_random.random() < params.re_haart_prob
            ):
                haart_params = (
                    self.agent.location.params.demographics[self.agent.race]
                    .sex_type[self.agent.sex_type]
                    .drug_type[self.agent.drug_type]
                    .haart
                )
                self.agent.haart.initiate(model.run_random, haart_params, "prob")
                self.re_initiated = True
  
        # stop tracing of this agent if time
        if self.active and model.time >= self.time + params.trace_duration:
            self.active = False
            self.time = None
        
        # this is to ensure agents are able to test negative multiple times throughout the simulation if contacted. On the otherhand, agents cannot test positive multiple times
        if self.tested_negative and model.time > self.tested_negative_time:
            self.tested_negative = False
            self.tested_negative_time = None
        
        # this is to ensure agents are able to participate multiple times throughout the simulation 
        if self.ps_participant and model.time > self.ps_participant_time:
            self.ps_participant = False
            self.ps_participant_time = None
    
    def trace_partners(self, params, model: "model.TITAN"):
        for ptnr in self.agent.get_partners(params.bond_type):
            
            # step 3
            self.info = True if model.run_random.random() < params.tracing_prob else False # can only trace if participant provides information
            # step 4
            self.contact = True if model.run_random.random() < params.contact_prob[ptnr.race] else False # contact stochastically
            
            if self.info:
                self.info_stat += 1
            
            if (
                self.info
                and self.contact
            ):
                self.contact_stat += 1
                ptnr.partner_tracing.active = True  # type: ignore[attr-defined]
                ptnr.partner_tracing.time = model.time  # type: ignore[attr-defined]
                
    def test_agent(self, agent_exposure, params, model: "model.TITAN"):
        
        # if agent has hiv, mark as diagnosed
        if (
            agent_exposure.active
            and model.run_random.random() < params.dx_prob
        ):    
            agent_exposure.diagnose(model) # markes agent with positive hiv diagnoses (not stochastic)
            self.ps_dx = True
            self.ps_dx_time = model.time
            self.add_positive_count()
            
        # if agent does not have hiv, mark as tested negative
        else:
            self.tested_negative = True
            self.tested_negative_time = model.time
            self.add_negative_count()
        
    def set_stats(self, stats: Dict[str, int], time: int):
        
        if self.participate:
            if self.participate_time == time:
                stats["step_1"] += 1
            
                stats["step_3"] += self.info_stat
                stats["step_4"] += self.contact_stat
        
        if self.ps_participant:
            if self.ps_participant_time == time:
                stats["ps_participants"] += 1
                
        if self.ps_dx:
            if self.ps_dx_time == time:
                stats["ps_dx"] += 1
            
        if self.tested_negative:
            if self.tested_negative_time == time:
                stats["ps_negative"] += 1
        
        if self.ps_prev_dx:
            if self.ps_prev_dx_time == time:
                stats["ps_prev_dx"] += 1
        if self.prev_dx_not_on_haart:
            stats["prev_dx_not_on_haart"] += 1
                    
        if self.re_initiated:
            stats["re_haart"] += 1
            
        stats["ps_negative_count"] = self.get_negative_count()
        
        stats["ps_positive_count"] = self.get_positive_count()
        