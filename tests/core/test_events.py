"""`linc/core/events.py` 的单元测试。"""
from datetime import datetime

from linc.core.events import ClusterCreated, DomainEvent, EventBus, RepresentativeAdded


def test_event_bus_publishes_to_all_subscribers_in_order():
    received = []
    bus = EventBus()
    bus.subscribe(lambda e: received.append(("sub1", e)))
    bus.subscribe(lambda e: received.append(("sub2", e)))

    event = ClusterCreated(timestamp=datetime.now(), payload={"cluster_id": "c1"})
    bus.publish(event)

    assert len(received) == 2
    assert received[0][0] == "sub1"
    assert received[1][0] == "sub2"
    assert received[0][1] is event


def test_event_bus_with_no_subscribers_does_not_raise():
    bus = EventBus()
    bus.publish(DomainEvent(timestamp=datetime.now(), payload={}))


def test_domain_event_payload_default_is_independent_dict():
    e1 = RepresentativeAdded(timestamp=datetime.now())
    e2 = RepresentativeAdded(timestamp=datetime.now())
    e1.payload["gene_id"] = "g1|1"
    assert "gene_id" not in e2.payload
